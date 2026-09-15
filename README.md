# **완전한 원형 Ring Resonator 기반 2채널 WDM ADD-Drop Filter 구현 및 시뮬레이션**  
  
# **1. 개요**  
이번에 알아볼 것은 이전 Meep 기반 원형/Racetrack Ring Resonator 분석에 이어서, **Tidy3D(3D FDTD)** 로 완전한 원형 Ring Resonator를 이용한 2채널 WDM Add-Drop 필터를 설계·검증하는 과정입니다.  
Ring Resonator의 스펙은 Bus와 Ring의 gap은 0.15μm이며 Ring의 radius는 5.0μm입니다. length_X, Y는 0이기에 이전 프로젝트와 마찬가지로 완전한 원형 Ring입니다. λ0 = 1.55μm이고 관심 대역은 1.50 ~ 1.60μm입니다.
waveguide 물질은 Si(n=3.47), Cladding 물질은 SiO2(n=1.44)이며, 시뮬레이션 도메인은 15×15×1.98μm입니다.

이번 프로젝트의 핵심 차이점은, 이전 원형 Ring Resonator에서 "coupling이 너무 약해 through만으로는 공진을 확인할 수 없었던" 문제를 gap을 0.2μm에서 0.15μm로 줄이고, 3D FDTD(Tidy3D)로 전환하여 재검증했다는 점입니다. 그리고 여기서 더 나아가, Ring이 하나의 through-port 소자가 아니라 **through 채널(비공진)과 drop 채널(공진) 두 개의 파장을 동시에 분리하는 WDM 필터**로 동작할 수 있는지를 확인하는 것이 목표입니다.

# **2. 이론적 예측 (시뮬레이션 이전)**

시뮬레이션을 돌리기 전에, 먼저 이 스펙에서 어떤 스펙트럼이 나올지, FSR은 몇일지를 이론적으로 먼저 예측해보았습니다.

## 2-1. FSR 예측

이전 프로젝트에서는 ng를 문헌값 근사(≈4.2)로 사용했지만, 이번에는 버스 도파로 단면(400×220nm strip)을 Tidy3D `ModeSolver`로 직접 풀어 neff(λ), ng(λ0)를 구했습니다.

$$n_g = n_{eff} - \lambda \frac{dn_{eff}}{d\lambda}$$

$$\text{FSR}(\lambda_0) = \frac{\lambda_0^2}{n_g \cdot L}, \qquad L = 2\pi R$$

| 항목 | 값 |
|---|---|
| neff(λ0=1.55μm) | 2.2168 |
| ng(λ0=1.55μm) | 4.2319 |
| L (둘레) | 31.4159 μm |
| **예측 FSR** | **18.071 nm** |

이는 이전 프로젝트에서 문헌 근사치(ng≈4.2)로 계산했던 18.2nm와 거의 비슷한 값이지만, 이번에는 근사가 아니라 실제 스펙(radius=5.0, width=0.4, thickness=0.22)에 대한 모드솔빙 결과이므로 더 신뢰할 수 있는 예측값입니다.

## 2-2. Add-Drop 스펙트럼 형태 예측 (TMM)

이전 프로젝트에서 사용했던 all-pass 링 공식

$$T = \frac{a^2 - 2ra\cos\theta + r^2}{1 - 2ra\cos\theta + (ra)^2}$$

은 through 포트 하나만 있는 all-pass 구조에 대한 식이었습니다. 이번에는 drop 포트가 추가된 add-drop 구조이므로, 버스-링 결합계수를 through측 t1, drop측 t2로 분리한 아래 형태를 사용합니다.

$$T_{through}=\left|\frac{t_1-t_2 a e^{j\phi}}{1-t_1 t_2 a e^{j\phi}}\right|^2,\qquad T_{drop}=\left|\frac{-\kappa_1\kappa_2\sqrt{a}\,e^{j\phi/2}}{1-t_1 t_2 a e^{j\phi}}\right|^2 \quad(\kappa_i=\sqrt{1-t_i^2})$$

t1, t2(결합), a(loss)는 gap=0.15μm만으로 해석적으로 바로 구하기는 어려우므로 별도의 directional coupler 모드해석이 필요합니다. 따라서 이 값들은 3장의 FDTD 결과를 얻은 뒤, 4장에서 역산하여 채워 넣었습니다. 즉 이번 프로젝트의 이론적 예측은 FSR은 시뮬레이션 전에 예측, 스펙트럼 형태는 FDTD 결과와 함께 검증하는 2단계로 진행했습니다.

# **3. 시뮬레이션 세팅**

Source는 `ModeSource`를 사용했고 λ0=1.55μm 중심으로 wavelength_width=0.1μm(1.50~1.60μm)의 GaussianPulse입니다. Through와 Drop 두 포트에 각각 `FluxMonitor`를 배치했고, 링 내부에 `FieldTimeMonitor`를 두어 정상상태 필드도 함께 관찰했습니다.

이전 Meep 기반 시뮬레이션과 다르게, 이번에는 별도의 reference straight waveguide의 정규화 과정 없이 Tidy3D의 flux 결과를 바로 dB로 환산해서 사용했습니다(`Decibels = 10·log10(P/P_in)`). grid는 `auto(wavelength=1.55, min_steps_per_wvl=15)`, `run_time=30ps`, 경계조건은 전면 PML입니다.

```python
modesource = td.ModeSource(
    center=(-5, -wg_center_y, 0), size=(0, 1.6, 1.32),
    source_time=td.GaussianPulse(freq0=freq0, fwidth=fwidth),
    direction='+', mode_spec=mode_spec, mode_index=0,
)
ThroughMonitor = td.FluxMonitor(center=(5, -wg_center_y, 0), size=(0, 1.6, 1.32), freqs=freqs, name='ThroughMonitor')
DropMonitor    = td.FluxMonitor(center=(-5, wg_center_y, 0), size=(0, 1.6, 1.32), freqs=freqs, name='DropMonitor')
```


# **4. 스펙트럼 분석 및 채널 선정**

Through/Drop 스펙트럼을 추출한 결과는 다음과 같습니다.

![FDTD Through 스펙트럼](./assets/fdtd_through.png)
![FDTD Drop 스펙트럼](./assets/fdtd_drop.png)

이전 원형 Ring Resonator(gap=0.2μm)에서는 through의 dip이 -0.02~-0.1dB 수준으로 거의 관측이 불가능한 수준이었는데, 이번 스펙(gap=0.15μm, 3D FDTD)에서는 FSR 간격으로 뚜렷한 dip(through 기준 최대 -21.8dB)이 여러 개 관측되었습니다. gap을 0.05μm 줄인 것만으로 결합 세기가 체감상 크게 달라졌고, 이전 결론("coupling length/gap이 결합 세기를 결정한다")과 같은 방향으로 재확인된 셈입니다.

이 스펙트럼에서 WDM 채널로 사용할 두 파장을 다음과 같이 선정했습니다.

| 채널 | 파장 | Through | Drop |
|---|---|---|---|
| λ1 (공진, drop 채널) | 1.56774 μm | −20.69708 dB | −0.2974766 dB |
| λ2 (비공진, through 채널) | 1.54048 μm | −0.0234321 dB | −24.55664 dB |

λ1은 through에서 깊게 꺼지고 drop으로 대부분 에너지가 전달되는 공진 파장, λ2는 반대로 through를 그대로 통과하고 drop으로는 거의 새어나가지 않는 비공진 파장입니다. 두 채널의 through/drop이 서로 상보적으로 동작하므로 2채널 WDM 필터로서 동작 조건을 만족합니다.

`ResonanceFinder`로 추출한 Q값은 다음과 같습니다.

| freq (Hz) | wavelength (μm) | Q |
|---|---|---|
| 1.845125e+14 | 1.624781 | 1077.0 |
| 1.867731e+14 | 1.605116 | 771.7 |
| 1.867738e+14 | 1.605109 | 1281.0 |
| 1.889927e+14 | 1.586265 | 1501.8 |
| **1.912390e+14** | **1.567633 (≈λ1)** | **1754.8** |
| 1.934702e+14 | 1.549554 | 2122.0 |
| 1.956994e+14 | 1.531903 | 2477.3 |
| 1.979503e+14 | 1.514484 | 2972.6 |
| 2.001742e+14 | 1.497658 | 3246.2 |
| 2.001827e+14 | 1.497594 | 3083.3 |
| 2.024382e+14 | 1.480909 | 4192.1 |

이전 원형 Ring Resonator의 Harminv 결과(Q=3112.8, 단일 모드)와 비교하면, 이번에는 대역 전체에서 여러 개의 모드가 규칙적으로 검출되었고 Q값도 1000~4000대로 합리적인 범위입니다. Q가 파장이 짧아질수록(고주파일수록) 커지는 경향도 보이는데, 이는 손실/결합 조건이 파장에 따라 서서히 변하기 때문으로 해석됩니다.

정상상태(t=3.22ps) 필드 분포도 확인했습니다. λ1, λ2를 동시에 여기했을 때 링 내부에 정상적인 공진 모드 형태가 형성됨을 직접 눈으로 확인할 수 있었습니다.

![정상상태 필드 분포](./assets/steady_state_field.png)

# **5. 예측값과 비교 검증**

## 5-1. FSR 검증

| | 값 |
|---|---|
| 2절 이론 예측 FSR | 18.071 nm |
| FDTD 실측 FSR (인접 공진 λ 간격, 1.567633→1.549554μm) | 18.079 nm |
| 오차 | **0.008 nm (≈0.04%)** |

Full 3D FDTD를 돌리기 전, 단면 모드솔빙만으로 예측한 FSR이 실제 FDTD 결과와 사실상 일치했습니다. 이전 Meep 프로젝트에서는 ng 근사치 사용으로 인한 오차가 한계점으로 지적되었는데, 이번에는 실제 스펙에 대한 모드솔빙을 사용해 그 한계를 개선한 결과로 볼 수 있습니다.

## 5-2. 독립 채널(λ2) 교차 검증

2-2절의 t1, t2, a는 λ1(공진점) 한 지점의 through/drop/Q 값으로만 역산(least-squares fitting)했습니다.

```python
sol = least_squares(residuals, x0=[0.9, 0.9, 0.999], bounds=(0,1))
# t1 = 0.93348, t2 = 0.92493, a = 0.995202
```

피팅에 사용하지 않은 λ2에서 이 모델이 예측한 값과 FDTD 실측값을 비교하면 다음과 같습니다.

| | Through | Drop |
|---|---|---|
| TMM 예측 (λ2) | −0.0251 dB | −22.6487 dB |
| FDTD 실측 (λ2) | −0.0234321 dB | −24.55664 dB |
| 차이 | 0.0017 dB | 1.91 dB |

Through는 사실상 일치합니다. Drop은 dB 기준 약 1.9dB 차이가 나는데, 왜 일까요? 선형 파워로 환산해보면 FDTD 0.00350, TMM 0.00543로 절대 파워 차이는 약 0.002(<0.3%) 수준입니다. 즉 이 지점은 drop이 -20dB 이하로 깊게 떨어지는 딥(dip) 근처이고, 이런 구간은 선형 파워의 작은 차이도 dB 스케일에서는 크게 증폭되어 보입니다. 이는 이전 프로젝트에서 "눈금(-15,1) 때문에 결과가 flat해 보였다"는 경험과 비슷한 맥락으로, 절대적인 물리량 차이보다 스케일 표현 방식이 차이를 과장할 수 있음을 다시 한번 보여주는 사례입니다.

# **6. 결론**
- 이전 원형 Ring Resonator(Meep, gap=0.2μm)에서는 coupling 구간이 짧아 through만으로 공진 여부를 확인하기 어려웠던 반면, 이번 스펙(Tidy3D, gap=0.15μm)에서는 through/drop 모두에서 뚜렷한 다중 공진 딥/피크가 확인되어 2채널 WDM Add-Drop 필터로 동작함을 확인했습니다.
- 시뮬레이션 이전에 단면 모드솔빙만으로 예측한 FSR(18.071nm)이 FDTD 실측(18.079nm)과 0.04% 오차로 일치했습니다.
- FDTD 공진점(λ1) 하나로 역산한 TMM 모델이, 피팅에 쓰이지 않은 독립 채널(λ2)에서도 선형 파워 기준 0.3% 이내로 일치하여 모델의 일반화 성능을 검증했습니다.
- 최종적으로 λ1=1.56774μm(drop 채널), λ2=1.54048μm(through 채널) 2채널 WDM 필터 설계를 완료했습니다.

## 한계점 및 아쉬운 점
- t1, t2, a는 gap으로부터 해석적으로 직접 유도한 것이 아니라 FDTD 결과를 역이용해 얻은 값입니다. gap→coupling 관계를 순수하게 예측하려면 별도의 directional coupler 모드해석(even/odd supermode) 단계가 필요합니다.
- λ1 캘리브레이션에 사용된 공진 파장이 이산 주파수 그리드상 정확한 공진점과 미세하게 어긋나 있을 수 있어, drop 채널처럼 dB 스케일 민감도가 큰 구간에서는 오차가 확대되어 보일 수 있습니다.

## 향후 방향
- gap을 변수로 하는 directional coupler 모드해석을 추가해 t1, t2를 gap으로부터 직접 예측하는 완전한 a priori 모델로 확장할 예정입니다.
- λ1 하나가 아니라 Q-table의 여러 공진 모드와 두 채널의 through/drop 값을 동시에 피팅하여 전 대역에서 더 안정적인 TMM 파라미터를 얻는 방향으로 개선할 계획입니다.

# 파일 구조

```
.
├── README.md
├── sim/
│   ├── ring_wdm_fdtd.py       # Tidy3D FDTD 본 시뮬레이션
│   └── mode_solve_theory.py   # neff/ng 모드솔빙 + TMM 이론 예측/피팅
├── data/
│   └── sim_data_ring.hdf5
└── assets/
    ├── fdtd_through.png
    ├── fdtd_drop.png
    └── steady_state_field.png
```


