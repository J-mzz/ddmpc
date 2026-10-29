# Robust Data-Driven Receding Horizon Control

MATLAB code for the paper

> J. Zheng, S. Kiani, M. Sznaier, and C. Lagoa, "Robust data-driven receding horizon control," *IFAC-PapersOnLine*, vol. 59, no. 16, pp. 25–30, 2025.
> [doi:10.1016/j.ifacol.2025.10.074](https://doi.org/10.1016/j.ifacol.2025.10.074) · [arXiv:2510.06153](https://arxiv.org/abs/2510.06153)

The controller handles an unknown linear system with bounded process noise. It describes all models consistent with noisy data by a polytope (the consistency set) and computes the largest robust invariant set inside the state constraints. It then applies a receding-horizon controller that refines the consistency set with every new measurement.

## Requirements

- MATLAB (tested with R2024a) with the Symbolic Math Toolbox and the Optimization Toolbox
- [YALMIP](https://yalmip.github.io/)
- [MOSEK](https://www.mosek.com/) (free academic licenses are available); tested with MOSEK 9.1

## Usage

Run either script from MATLAB. Each one adds `utils/` to the path, generates the data, computes the invariant set, simulates the closed loop with and without receding horizon, and plots the results.

```matlab
example1   % Example 1 (Fig. 1): noise bound 0.04, random disturbance
example2   % Example 2 (Fig. 2): noise bound 0.1, constant disturbance
```

| Script | Paper | Settings | Figures |
|---|---|---|---|
| `example1.m` | Example 1, Fig. 1 | `eps = 0.04`, `x0 = [1; 0.8]`, 50 steps | (a) largest invariant set, (b) trajectories, (c) contraction rate λ, (d) control input |
| `example2.m` | Example 2, Fig. 2 | `eps = 0.1`, `x0 = [-1; 0.8]`, 20 steps | same as Fig. 1 |

The settings are at the top of each script. Both examples use the system `A = [0 -0.99; 0.99 0]`, `B = [0; 1]` (used to generate the data and to simulate the closed loop, see `utils/getSystem.m`), the state constraint `||x||_inf <= 1`, the input constraint `|u| <= 1`, and 10 initial noisy samples.

## Repository layout

```
example1.m, example2.m     scripts for Examples 1 and 2
utils/
  ddmpc.m                  one step of the data-driven controller (an LP)
  maxis.m                  one iteration of the largest invariant set computation
  getDataRND.m             noisy samples of the true system
  getCons.m                data matrices of the consistency set
  updateDyn*.m             one simulation step (random or constant disturbance)
  getSystem.m, setCons.m   system definition and problem dimensions
  third_party/             external code, see utils/third_party/README.md
```

## Citation

If you use this code, please cite:

```bibtex
@article{zheng2025robust,
  title={Robust data-driven receding horizon control},
  author={Zheng, Jian and Kiani, Sahand and Sznaier, Mario and Lagoa, Constantino},
  journal={IFAC-PapersOnLine},
  volume={59},
  number={16},
  pages={25--30},
  year={2025},
  publisher={Elsevier},
  doi={10.1016/j.ifacol.2025.10.074}
}
```

## License

The code in this repository is released under the [MIT License](LICENSE), except the files in `utils/third_party/`, which keep their original licenses.
