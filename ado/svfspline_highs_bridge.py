#!/usr/bin/env python3
"""Optional HiGHS backend for the fixed-hyperparameter SVF-Spline QP.

This bridge is intentionally separate from the default Mata implementation.
It accepts only CSV files supplied by svfspline.ado and writes output files
only after HiGHS reports an optimal solution and independently recomputed
primal and KKT diagnostics meet the requested tolerances.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import cvxpy as cp
import numpy as np


def read_csv(path: Path) -> np.ndarray:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
    if data.dtype.names is None:
        raise ValueError(f"CSV has no header: {path}")
    return np.column_stack([data[name] for name in data.dtype.names])


def design(x: np.ndarray, d: int, ranges: np.ndarray) -> np.ndarray:
    blocks = []
    for j in range(x.shape[1]):
        lo, hi = ranges[j]
        nodes = np.linspace(lo, hi, d + 1)
        blocks.append(np.column_stack((np.ones(x.shape[0]), x[:, j],
                                      *[np.maximum(x[:, j] - knot, 0.0)
                                        for knot in nodes])))
    return np.column_stack(blocks)


def constraints(phi: np.ndarray, y: np.ndarray, d: int, eps: float, n_inputs: int,
                beta: cp.Variable, xi: cp.Variable):
    fitted = phi @ beta
    cons = [y - fitted <= 0, fitted - xi - y - eps <= 0, -xi <= 0]
    q = phi.shape[1]
    for j in range(n_inputs):
        off = j * (d + 3)
        unit = np.zeros(q)
        unit[off] = -1
        cons.append(unit @ beta <= 0)
        unit = np.zeros(q)
        unit[off + 1] = -1
        cons.append(unit @ beta <= 0)
        for k in range(off + 2, off + d + 3):
            unit = np.zeros(q)
            unit[k] = 1
            cons.append(unit @ beta <= 0)
        for k in range(2, d + 4):
            unit = np.zeros(q)
            unit[(off + 1):(off + k)] = -1
            cons.append(unit @ beta <= 0)
    return fitted, cons


def solve_one(phi: np.ndarray, y: np.ndarray, d: int, c_pen: float, eps: float):
    n, q = phi.shape
    beta = cp.Variable(q)
    xi = cp.Variable(n)
    fitted, cons = constraints(phi, y, d, eps, q // (d + 3), beta, xi)
    problem = cp.Problem(cp.Minimize(cp.sum_squares(beta) + c_pen * cp.sum(xi)), cons)
    value = problem.solve(solver=cp.HIGHS)
    if problem.status != cp.OPTIMAL:
        raise RuntimeError(f"HiGHS status: {problem.status}")

    fitted_value = np.asarray(fitted.value).reshape(-1)
    theta = np.r_[np.asarray(beta.value).reshape(-1), np.asarray(xi.value).reshape(-1)]
    if not np.all(np.isfinite(theta)):
        raise RuntimeError("HiGHS returned non-finite primal values")
    if not np.all(np.isfinite(fitted_value)):
        raise RuntimeError("HiGHS returned non-finite fitted values")
    # Reconstruct A theta <= b in the same order as constraints().
    rows, rhs, dual = [], [], []
    rows.extend(np.c_[-phi, np.zeros((n, n))]); rhs.extend(-y); dual.extend(cons[0].dual_value)
    rows.extend(np.c_[phi, -np.eye(n)]); rhs.extend(y + eps); dual.extend(cons[1].dual_value)
    rows.extend(np.c_[np.zeros((n, q)), -np.eye(n)]); rhs.extend(np.zeros(n)); dual.extend(cons[2].dual_value)
    for con in cons[3:]:
        # Shape restrictions contain beta only; CVXPY exposes their coefficient
        # via the construction below, rebuilt from d rather than introspection.
        dual.extend(np.asarray(con.dual_value).reshape(-1))
    shape_rows = []
    for j in range(q // (d + 3)):
        off = j * (d + 3)
        for index, sign in ((off, -1), (off + 1, -1)):
            row = np.zeros(q + n); row[index] = sign; shape_rows.append(row)
        for index in range(off + 2, off + d + 3):
            row = np.zeros(q + n); row[index] = 1; shape_rows.append(row)
        for k in range(2, d + 4):
            row = np.zeros(q + n); row[(off + 1):(off + k)] = -1; shape_rows.append(row)
    rows.extend(shape_rows); rhs.extend(np.zeros(len(shape_rows)))
    a = np.vstack(rows)
    b = np.asarray(rhs)
    lam = np.asarray(dual, dtype=float)
    if len(lam) != len(b):
        raise RuntimeError("HiGHS dual-vector dimension mismatch")
    slack = b - a @ theta
    gradient = np.r_[2 * theta[:q], np.full(n, c_pen)]
    stationarity = gradient + a.T @ lam
    scale = max(1.0, np.max(np.abs(gradient)), np.max(np.abs(a.T @ lam)))
    diagnostics = {
        "objective": float(value),
        "max_violation": float(max(0.0, -np.min(slack))),
        "max_stationarity": float(np.max(np.abs(stationarity)) / scale),
        "max_complementarity": float(np.max(np.abs(lam * slack)) / max(1.0, abs(value))),
        "iterations": int(problem.solver_stats.num_iters or 0),
    }
    if (not all(np.isfinite(value) for value in diagnostics.values()) or
            diagnostics["max_violation"] > 1e-7 or diagnostics["max_stationarity"] > 1e-7 or
            diagnostics["max_complementarity"] > 1e-7):
        raise RuntimeError(f"HiGHS certification failed: {diagnostics}")
    return theta[:q], fitted_value, diagnostics


def cross_validate(x: np.ndarray, y: np.ndarray, cgrid: list[float],
                   epsgrid: list[float], dgrid: list[int], folds: int,
                   seed: int):
    """Paper-aligned, leakage-free CV with independently certified QP fits."""
    n, n_outputs = y.shape
    rng = np.random.default_rng(seed)
    fold_id = np.empty(n, dtype=int)
    fold_id[rng.permutation(n)] = np.arange(n) % folds + 1
    rows, fold_rows, pred_rows = [], [], []
    for d in dgrid:
        for c_pen in cgrid:
            for eps in epsgrid:
                losses, failed = [], False
                for fold in range(1, folds + 1):
                    train = fold_id != fold
                    valid = ~train
                    ranges = np.c_[x[train].min(axis=0), x[train].max(axis=0)]
                    phi_train = design(x[train], d, ranges)
                    phi_valid = design(x[valid], d, ranges)
                    squared_error = 0.0
                    fold_predictions = []
                    try:
                        for output in range(n_outputs):
                            beta, _, _ = solve_one(phi_train, y[train, output], d, c_pen, eps)
                            prediction = phi_valid @ beta
                            if not np.all(np.isfinite(prediction)):
                                raise RuntimeError("HiGHS CV prediction is non-finite")
                            squared_error += float(np.sum((y[valid, output] - prediction) ** 2))
                            for obs, value in zip(np.flatnonzero(valid), prediction):
                                fold_predictions.append((c_pen, eps, d, fold, obs + 1, output + 1, value))
                    except Exception:
                        failed = True
                        fold_rows.append((c_pen, eps, d, fold, np.nan))
                        break
                    loss = squared_error / (int(valid.sum()) * n_outputs)
                    losses.append(loss)
                    fold_rows.append((c_pen, eps, d, fold, loss))
                    pred_rows.extend(fold_predictions)
                rows.append((c_pen, eps, d, np.nan if failed else float(np.mean(losses))))
    cv = np.asarray(rows, dtype=float)
    eligible = np.flatnonzero(np.isfinite(cv[:, 3]))
    if not len(eligible):
        raise RuntimeError("all HiGHS CV candidates failed certification")
    best = eligible[0]
    for index in eligible[1:]:
        if cv[index, 3] < cv[best, 3] - 1e-14:
            best = index
    return cv, np.asarray(fold_rows, dtype=float), np.asarray(pred_rows, dtype=float), fold_id, cv[best]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--x", type=Path, required=True)
    p.add_argument("--y", type=Path, required=True)
    p.add_argument("--beta", type=Path, required=True)
    p.add_argument("--frontier", type=Path, required=True)
    p.add_argument("--ranges", type=Path, required=True)
    p.add_argument("--diagnostics", type=Path, required=True)
    p.add_argument("--d", type=int, required=True)
    p.add_argument("--c", type=float, required=True)
    p.add_argument("--eps", type=float, required=True)
    p.add_argument("--cgrid")
    p.add_argument("--epsgrid")
    p.add_argument("--dgrid")
    p.add_argument("--folds", type=int)
    p.add_argument("--seed", type=int)
    p.add_argument("--cvresults", type=Path)
    p.add_argument("--cvfold", type=Path)
    p.add_argument("--cvpred", type=Path)
    p.add_argument("--foldid", type=Path)
    a = p.parse_args()
    x, y = read_csv(a.x), read_csv(a.y)
    if x.shape[0] != y.shape[0]:
        raise ValueError("X and Y have different row counts")
    if any(value is not None for value in (a.cgrid, a.epsgrid, a.dgrid, a.folds, a.seed,
                                            a.cvresults, a.cvfold, a.cvpred, a.foldid)):
        if not all(value is not None for value in (a.cgrid, a.epsgrid, a.dgrid, a.folds, a.seed,
                                                   a.cvresults, a.cvfold, a.cvpred, a.foldid)):
            raise ValueError("HiGHS CV arguments must be supplied together")
        cv, cvfold, cvpred, foldid, selected = cross_validate(
            x, y, [float(v) for v in a.cgrid.split()],
            [float(v) for v in a.epsgrid.split()], [int(v) for v in a.dgrid.split()],
            a.folds, a.seed)
        a.c, a.eps, a.d = float(selected[0]), float(selected[1]), int(selected[2])
        np.savetxt(a.cvresults, cv, delimiter=",", header="C,eps,d,mse", comments="")
        np.savetxt(a.cvfold, cvfold, delimiter=",", header="C,eps,d,fold,mse", comments="")
        np.savetxt(a.cvpred, cvpred, delimiter=",", header="C,eps,d,fold,obs,output,prediction", comments="")
        np.savetxt(a.foldid, foldid[:, None], delimiter=",", header="fold", comments="")
    ranges = np.c_[x.min(axis=0), x.max(axis=0)]
    phi = design(x, a.d, ranges)
    beta, frontier = [], []
    diagnostics = []
    for output in range(y.shape[1]):
        b, f, diag = solve_one(phi, y[:, output], a.d, a.c, a.eps)
        beta.append(b); frontier.append(f); diagnostics.append(diag)
    beta_arr, frontier_arr = np.column_stack(beta), np.column_stack(frontier)
    summary = {
        "objective": sum(d["objective"] for d in diagnostics),
        "max_violation": max(d["max_violation"] for d in diagnostics),
        "max_stationarity": max(d["max_stationarity"] for d in diagnostics),
        "max_complementarity": max(d["max_complementarity"] for d in diagnostics),
        "iterations": max(d["iterations"] for d in diagnostics),
    }
    for path, data, prefix in ((a.beta, beta_arr, "b"), (a.frontier, frontier_arr, "f"),
                               (a.ranges, ranges, "r")):
        np.savetxt(path, data, delimiter=",", header=",".join(f"{prefix}{i+1}" for i in range(data.shape[1])), comments="")
    np.savetxt(a.diagnostics, np.array([[summary[k] for k in summary]]), delimiter=",",
               header=",".join(summary), comments="")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"SVFSPLINE_HIGHS_ERROR: {exc}", file=sys.stderr)
        raise
