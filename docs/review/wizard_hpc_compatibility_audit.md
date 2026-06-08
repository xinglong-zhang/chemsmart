# wizard HPC OS/scheduler compatibility audit

검토 기준 커밋: `f58eb046` (`fix(agent): Python 3.10 compat + SGE wizard auto-discovery + provider UX (#81)`)  
검토 방식: 코드 리뷰 only, 소스 수정 없음.  
읽은 파일: `chemsmart/agent/wizard/scheduler_env.py`, `probe.py`, `parsers.py`, `survey.py`, `normalize.py`, `software.py`, `providers.py`, `docs/installation-troubleshooting.md`.  
확인 실행: `pytest -v tests/agent/wizard/test_parsers.py tests/agent/wizard/test_survey.py tests/agent/wizard/test_probe.py tests/agent/wizard/test_software.py` → `23 passed`.

## 1. Verdict

- **가장 가능성이 높은 성공 조합**은 **SoGE 8.1 / 일부 UGE 8.6+ / SLURM 20-23** + **EL8/EL9 또는 Ubuntu 22.04** + **bash 로그인 셸** + **scheduler 바이너리가 이미 `/usr/bin` 또는 PR #81 후보 경로 안에 있는 경우**다 (`scheduler_env.py:28-57`, `probe.py:227-259`).
- **실운영에서 바로 깨질 가능성이 높은 조합**은 **LSF 9/10**, **Torque 6**, **Bright/xCAT 계열의 `/cm/local/apps/*/current` 설치**, **tcsh 기본 셸**, **Spack-only 소프트웨어 스택**이다. 지금 코드는 여기를 거의 안 본다.
- **SGE 계열은 “1번 클러스터”에서는 통과했지만**, 그 성공은 `qhost` fallback이 먹히는 **SoGE 8.1.8 / CentOS 7 / local mode A**에 강하게 편향돼 있다. **queue-wide `slots`, alarm 상태, mixed-node queue**에서는 조용히 틀릴 수 있다 (`parsers.py:235-306`, `survey.py:159-220`).
- **fix #1, #2, #8은 대체로 OS 중립적**이고, 이번 확장 감사의 핵심 위험은 **fix #5-#7**(scheduler auto-discovery, qhost fallback, slot 보정)의 범위 한계다.
- **Cray PALS / ALPS는 사실상 미지원**이다. 운 좋게 underlying SLURM/PBS가 표준 PATH에 있으면 일부가 우연히 동작할 수 있지만, 그건 PALS/ALPS 지원이 아니다.

## 2. Compatibility matrix

| Scheduler | auto-discovery path coverage | parser correctness | bin layout assumption | gaps |
|---|---|---|---|---|
| SGE 6.2u5 | ⚠️ `SGE_ROOT` 또는 `/opt/sge`, `/usr/sge`, `/usr/local/sge`면 가능하지만 `/cm/local/apps/sge/current` 류는 누락 (`scheduler_env.py:28-39`, `81-106`). | ⚠️ `qconf -sql/-sq`와 `qhost` 형식은 대체로 맞지만 `slots→cores`, state, mixed-node 처리 가정이 약하다 (`parsers.py:235-306`). | ⚠️ `util/arch` + `<cell>/common/settings.sh` 전제를 둔다; wrapper install이나 비표준 cell naming이면 빗나간다 (`scheduler_env.py:64-139`). | queue-wide `slots`, alarm/suspend 상태, remote mode B 무보강. |
| OGS 2011.11 | ⚠️ SGE와 동일 계열이어서 기본 루트면 가능하나 `current` symlink 계열 누락. | ⚠️ 출력 형식은 유사하나 오래된 OGS에서 queue state/complex naming 편차를 못 흡수한다. | ⚠️ `settings.sh`/`util/arch` 구조를 그대로 기대한다. | fix #6/#7이 serial queue에는 도움되지만 오래된 hostgroup/PE 관례는 여전히 공백. |
| SoGE 8.1 | ✅ PR #81 검증 대상이 바로 SoGE 8.1.8/CentOS7이며 local mode A에서는 가장 유리하다. | ⚠️ `h_vmem=INFINITY`/`slots=1`은 잡았지만 그 외 `slots>1`, `state=alarm`은 그대로다 (`docs/installation-troubleshooting.md:160-252`). | ⚠️ `/opt/sge`, `/cm/shared/apps/sge` 류는 맞지만 Bright의 `/cm/local/apps/sge/current`는 빠진다. | “1번 클러스터” 편향이 가장 큰 행. |
| UGE 8.6+ | ⚠️ `/opt/uge`, `/opt/univa`, `/cm/shared/apps/uge`는 보지만 `/cm/local/apps/uge/current`, `/cm/shared/apps/uge/current`는 놓친다 (`scheduler_env.py:28-39`). | ⚠️ `qconf`/`qhost` 자체는 읽히겠지만 UGE에서 흔한 queue-wide `slots` 합계는 per-job core로 오인할 수 있다. | ⚠️ Univa/Altair 패키징의 `current` symlink, Bright export path를 가정하지 않는다. | UGE에서 더 흔한 large shared queue에서 `NUM_CORES` 과대기재 위험. |
| SLURM 18/20 | ⚠️ `/usr/bin`, `/opt/slurm/bin`, `/usr/local/slurm/bin`면 되나 `/opt/slurm/current/bin`, `/cm/local/apps/slurm/current/bin`, OHPC 경로는 빠진다 (`scheduler_env.py:41-48`). | ⚠️ `sinfo --json`이 실패해도 `scontrol show partition --oneliner` 폴백은 있다; 다만 `MaxCPUsPerNode`/`TotalCPUs`를 per-job core처럼 쓰는 건 거칠다 (`survey.py:80-105`, `parsers.py:66-104`). | ⚠️ “bin만 PATH에 있으면 충분”하다고 보고 daemon/env skew는 무시한다. | 버전 skew, old field names, shared-node partition 의미 차이 미반영. |
| SLURM 22/23 | ✅ distro package 또는 표준 `/usr/bin` 설치에서는 경로 발견률이 높다. | ⚠️ JSON/oneliner 둘 다 있을 가능성이 높아도 `DefMemPerCPU`, `TRES`, heterogeneous partition 해석은 여전히 단순하다. | ⚠️ `/cm/local/apps/slurm/current/bin`, Cray/ Bright current path는 여전히 미포함. | 기본 동작은 제일 낫지만 mixed partition resource semantics는 여전히 위험. |
| LSF 9/10 | ❌ `scheduler_env.py`가 LSF를 전혀 보강하지 않는다; `bqueues`가 PATH에 이미 있어야만 탐지된다. | ❌ `bqueues -l` parser가 `PROCLIMIT`를 core 수처럼 읽고 `RUNLIMIT` 분/초 표기를 거의 못 읽는다 (`parsers.py:193-232`). | ❌ 보편 경로인 `/opt/ibm/lsfsuite/lsf/<ver>/<arch>/bin` 및 `LSF_ENVDIR`/`LSF_SERVERDIR` 가정이 없다. | 실제 운영 LSF는 path/env/bootstrap, parser 둘 다 취약. |
| OpenPBS 19+ | ⚠️ `/opt/pbs/bin`은 보지만 흔한 `/opt/pbs/default/bin`, Bright의 `current` symlink는 빠진다 (`scheduler_env.py:50-57`). | ⚠️ JSON/text 둘 다 시도하는 점은 좋지만 `select=`/`default_chunk`/site resource 이름 다양성은 충분히 못 읽는다 (`parsers.py:107-190`). | ⚠️ bin만 찾으면 된다고 보고 `PBS_EXEC` 기반 path 변형을 충분히 반영하지 않는다. | 표준 workq는 될 수 있으나 chunk/select 중심 사이트는 facts 누락 가능. |
| PBS Pro 18+ | ⚠️ OpenPBS와 유사; `/opt/pbs/bin` 표준이면 가능. | ⚠️ 기본 큐 판단, `resources_default.*` 중심 파싱은 PBS Pro의 site policy 차이를 많이 흘린다. | ⚠️ `/opt/pbs/default/bin`, vendor wrapper path를 모른다. | 특히 GPU/resource group naming 차이를 못 읽는다. |
| Torque 6 | ⚠️ `/opt/torque/bin`, `/cm/shared/apps/torque/current/bin`은 보지만 나머지 legacy 경로는 제한적. | ❌ JSON probe는 무의미하고 text parser도 `nodes=1:ppn=16`, `resources_default.nodes` 관례를 이해하지 못한다. | ⚠️ Torque는 오래된 package layout 변종이 많아 현재 후보군이 얕다. | legacy PBS row 중 가장 취약. |
| Cray PALS | ❌ PALS 전용 probe/path가 없다; underlying SLURM가 우연히 표준 PATH에 있어야만 일부 동작. | ❌ parser도 PALS 자체를 이해하지 못하고 SLURM facts만 본다. | ❌ `/opt/cray/pe/pals/default/bin` 류를 전혀 고려하지 않는다. | “PALS 지원”이 아니라 “표준 SLURM가 보이면 일부 우연히 통과” 수준. |
| Cray ALPS | ❌ ALPS 전용 probe/path가 없다; PBS/Moab layer가 우연히 보이면 일부만 동작. | ❌ `aprun`, ALPS placement, Cray login wrapper 정보는 전혀 수집하지 않는다. | ❌ `/opt/cray/alps/default/bin` 류 전제 부재. | 사실상 unknown/unsupported. |

## 3. OS audit

| OS | 기본 Python / Py3.10 포함 여부 | GCC ≥ C99 관점 | PR #81 후보군과 어긋나는 흔한 scheduler 경로 | 기본 셸 / 로그인 프로파일 포인트 |
|---|---|---|---|---|
| RHEL/CentOS 7 | 기본은 Python 2.7, `python3`도 보통 3.6 미만; **Py3.10 기본 미포함**. fix #1의 이식성 이슈가 가장 현실적이다. | 기본 GCC 4.8 계열이라 C99는 옵션으로 가능해도 wheel이 없으면 fix #2 상황이 바로 재현된다. | `/cm/local/apps/uge/current`, `/cm/local/apps/slurm/current/bin`, `/opt/ibm/lsfsuite/lsf/.../bin`, `/opt/pbs/default/bin`. | bash가 흔하지만 site module init이 `.bashrc`나 `/etc/bashrc` 전용이면 `bash -lc`만으로는 안 잡힌다. |
| RHEL/CentOS 8 | 기본 `python3`는 3.6 계열, **Py3.10 기본 미포함**. | GCC 8 계열 패키지는 C99 문제를 크게 줄이지만, compiler 미설치 노드는 여전히 빌드 취약. | `/opt/slurm/current/bin`, `/opt/pbs/default/bin`, `/cm/local/apps/*/current`. | `/etc/profile.d/*.sh` 기반 init은 잘 맞지만 non-bash default shell 사이트는 동일하게 취약. |
| RHEL/CentOS 9 | 기본 `python3`는 3.9 계열, **Py3.10 기본 미포함**. | GCC 11 계열이라 fix #2 중요도는 낮아지지만 source build 회피 전략 자체는 여전히 유효. | `/usr/lib64/slurm`, `/opt/slurm/current/bin`, `/cm/local/apps/slurm/current/bin`. | bash login은 대체로 무난; zsh/tcsh 사용자 계정이면 module init 비대칭은 남는다. |
| Rocky/Alma 8 | RHEL8과 사실상 동일; **Py3.10 기본 미포함**. | GCC 8 계열. site repos/EasyBuild/Spack overlay를 더 많이 탄다. | Bright에서 특히 `/cm/local/apps/uge/current`, `/cm/local/apps/slurm/current/bin`가 흔하다. | Bright 이미지에서 `/etc/profile.d`가 깔끔하면 되지만, 사용자 dotfile 의존 사이트는 깨진다. |
| Rocky/Alma 9 | RHEL9과 동일; **Py3.10 기본 미포함**. | GCC 11 계열. | `/opt/pbs/default/bin`, `/cm/local/apps/pbspro/current/bin`, `/opt/ibm/lsfsuite/...`. | bash login은 무난하지만 JupyterHub/OOD에서 non-login shell이면 probe 결과가 달라질 수 있다. |
| Ubuntu 20.04 LTS | 기본 `python3`는 3.8, **Py3.10 기본 미포함**. | GCC 9 계열이라 C99 빌드 조건은 충분하다. | distro package면 `/usr/bin`, custom build면 `/opt/slurm/current/bin`, `/opt/pbs/default/bin`. | `/etc/profile` + `~/.profile` 경로가 깔끔하면 bash -lc가 잘 맞는다. |
| Ubuntu 22.04 LTS | 기본 `python3`는 3.10, **Py3.10 기본 포함**. PR #81의 Python 호환성 검증 OS로 가장 편하다. | GCC 11 계열. fix #2는 상대적으로 낮은 위험. | module-only cluster라도 scheduler는 `/usr/bin`에 설치되는 경우가 많아 SLURM 쪽 성공률이 높다. | Lmod/env-modules init을 `/etc/profile.d`에 넣었다면 잘 맞고, user alias `ml` 의존이면 안 맞는다. |
| Ubuntu 24.04 LTS | 기본 `python3`는 3.12, **Py3.10 아님**. | GCC 13 계열. build 툴체인은 가장 여유 있다. | custom vendor stack이면 `/opt/slurm/current/bin`, `/opt/pbs/default/bin`, `/opt/ibm/lsfsuite/...`. | bash login은 무난하지만 zsh 사용자 비중이 높아 `bash -lc` 비대칭이 더 눈에 띈다. |
| SLES 15 | 기본 Python은 대체로 3.6 계열(서비스팩/모듈에 따라 상향 가능)로 **Py3.10 baseline 아님**. | GCC 7+ SDK/compilers는 C99에 충분하나, base install만 믿으면 compiler 부재가 잦다. | `/opt/slurm/current/bin`, `/opt/pbs/default/bin`, `/opt/ibm/lsfsuite/...`, `/opt/cray/pe/lmod`, `/opt/cray/pe/pals/default/bin`. | **tcsh 기본 계정이 아직 많다.** `probe.py`의 `bash -lc` 고정은 여기서 특히 잘 드러난다. |
| Cray Linux / Cray service OS | login/service node는 SLES 파생이 많고 compute node는 별도 이미지다; **Py3.10 기본 전제 금지**. | Cray PE는 `/usr/bin/gcc`보다 `cc/CC/ftn` wrapper가 실사용 경로다. fix #2의 “system GCC” 가정은 그대로 옮기기 어렵다. | `/opt/cray/pe/pals/default/bin`, `/opt/cray/alps/default/bin`, `/opt/slurm/default/bin`, site wrapper path. | module init이 Cray PE/Lmod 스크립트에 강하게 묶여 있어, login shell 종류/launch context에 따라 결과가 크게 달라진다. |

보충 메모:

- `providers.py`의 fix #8 (`load_dotenv` 선행)는 pure-Python reorder라서 위 OS 차이의 핵심 리스크는 아니다 (`providers.py:157-193`).
- `docs/installation-troubleshooting.md:246-252`의 “수동 scheduler load는 보통 필요 없음” 문장은 **CentOS 7 + local mode + candidate path 안의 SGE**에는 맞지만, 위 표의 다수 조합에는 그대로 일반화하면 안 된다.

## 4. Module system audit

### 4.1 현재 코드가 실제로 하는 일

- module system detection entry point는 `detect_module_system()`이며, 순서는 `type module` → 실패 시 `which module` → 성공하면 `module --version`이다 (`software.py:42-67`).
- 이 probe들은 모두 `ProbeRunner`를 통해 실행되고, local/remote 모두 **강제로 `bash -lc`**로 감싼다 (`probe.py:216-267`).
- module candidate 검색은 `module -t avail` 출력에서 정규식 매칭 줄만 뽑아 첫 토큰을 취하는 정도다 (`software.py:202-233`).

### 4.2 env-modules(Tcl)

- 대부분의 env-modules는 `/usr/share/Modules/init/bash` 또는 `/etc/profile.d/modules.sh`가 **로그인 셸에서 shell function `module`** 을 정의해 줘야 한다.
- 따라서 `type module` (`probe.py:140-143`)은 **function init이 이미 끝난 bash login shell** 에서는 맞지만, 다음 조건 중 하나면 곧바로 false negative가 난다.
  - 사용자가 tcsh/zsh를 기본 셸로 쓰고 bash 프로파일에 module init이 없음.
  - site가 `module` init을 `.bashrc`에만 넣고 `.bash_profile`에서는 안 부름.
  - OOD/JupyterHub/IDE launch가 non-login shell로 시작해 parent env가 비어 있음.
- `which module` fallback (`probe.py:145-149`, `software.py:50-57`)도 믿기 어렵다. 외부 `which`는 shell function을 못 보고, bash alias 기반 `which`도 non-interactive에서 비활성일 수 있다.

### 4.3 Lmod

- Lmod도 실체는 shell function/alias bootstrap이므로 첫 실패 원인은 env-modules와 거의 같다.
- `module --version`은 Lmod에서 stderr로 버전을 뿌리는 경우가 흔한데, 이건 `_merge_output()`이 흡수하므로 괜찮다 (`software.py:59-67`, `264-269`).
- 하지만 현 코드는 `ml` alias를 전혀 보지 않고, `module -t avail` 결과에 붙는 `(D)`, hidden/default marker, modulepath header, localized noise를 최소한만 지운다 (`software.py:225-233`). 그래서 **“있는데 없는 것처럼 보이는”** 경우보다 **“있는데 후보 문자열이 지저분하게 들어오는”** 경우가 다음 단계 문제다.

### 4.4 Spack

- Spack는 여기서 사실상 **0점**이다.
- `spack load`, `spack env activate`, `spack find --format`, `spack location -i` 중 아무 것도 probe catalog에 없다 (`probe.py:50-208`).
- 결과적으로 module system이 없는 Spack-only 사이트에서는 `ModuleSystem(kind="none")` + `ProgramFinding(source="none")`로 떨어지며, wizard는 **소프트웨어가 없는 환경**으로 오판한다.

### 4.5 왜 `software.type_module = ("type", "module")`가 클러스터별로 실패하는가

1. `type` 자체는 bash builtin이라 `bash -lc` 안에서는 실행되지만, **`module`이 bash function으로 준비돼 있지 않으면 바로 실패**한다.
2. **“로그인 셸”과 “인터랙티브 셸”은 다르다.** `bash -lc`는 `/etc/profile`과 `~/.bash_profile`은 읽지만 `.bashrc`는 보장하지 않는다.
3. **기본 셸이 tcsh/zsh인 사용자**는 site가 bash init을 별도로 안 넣어 두면 `module`이 생기지 않는다.
4. **`which module`은 shell function 탐지기가 아니다.** env-modules/Lmod에서 가장 흔한 정의 형태와 어긋난 fallback이다.
5. **remote mode B는 scheduler_env 같은 보강이 전혀 없다.** 즉 local mode A보다 module false negative가 더 자주 난다 (`probe.py:227-234` vs `248-267`).
6. **Spack/EasyBuild wrapper-only 사이트는 아예 모델 밖**이다. “module을 못 찾음”이 곧 “software stack 없음”은 아니다.

### 4.6 운영적 해석

- 지금 구조에서 **scheduler 발견 성공**과 **software stack 발견 성공**은 별개다.
- PR #81은 scheduler path 쪽만 일부 자동화했고 (`scheduler_env.py:243-255`), module/software 쪽은 여전히 **셸 초기화 품질에 전적으로 의존**한다.
- 그래서 “wizard가 scheduler는 찾는데 Gaussian/ORCA module은 못 찾는” 2번 클러스터가 매우 쉽게 나온다.

## 5. Numbered gap list (priority = likelihood × impact)

### 1. SGE/UGE queue-wide `slots`를 per-job core처럼 사용
- 관련 fix: **#7 부분 보완 실패**. `slots <= 1`일 때만 `qhost` CPU로 대체한다 (`survey.py:201-205`).
- 트리거: UGE 8.6+, SoGE shared queue, PE queue, large hostgroup queue.
- 제안: `qconf -sq`의 `slots`를 바로 core로 쓰지 말고, `qhost -q` 또는 queue-instance/hostgroup mapping을 추가해 **queue별 min/max slot fact**를 따로 계산.
- 심각도: **P0 silently-wrong**.

### 2. SGE state parser가 `disabled`만 본다
- 관련 fix: **#6/#7과 무관하게 놓친 부분**. `enabled/started`가 `"disabled" not in state.lower()`뿐이다 (`parsers.py:291-306`).
- 트리거: alarm/suspend/unknown/subordinate 상태가 흔한 production SGE/UGE.
- 제안: `state` 토큰을 `d`,`D`,`a`,`A`,`s`,`u`,`E` 등으로 해석하고, 선택 로직에서 unusable queue를 제외.
- 심각도: **P0 silently-wrong**.

### 3. SLURM partition core 수가 `MaxCPUsPerNode`/`TotalCPUs`로 과대 추정될 수 있음
- 관련 fix: **#5가 새 survey path를 열었지만 의미 해석은 단순함** (`parsers.py:78-87`).
- 트리거: shared-node partition, oversubscribe, heterogeneous partition, `DefCpuPerGPU`/QoS 제한 중심 사이트.
- 제안: `DefaultCPUsPerNode` 우선, 없으면 `DefMemPerCPU`와 함께 conservative fallback; `TotalCPUs`는 마지막 최후 수단으로만 사용.
- 심각도: **P0 silently-wrong**.

### 4. PBS/OpenPBS/PBS Pro/Torque에서 `select`/`nodes:ppn` 계열 자원을 못 읽음
- 관련 fix: **#5 범위 밖으로 남음** (`parsers.py:113-143`, `147-190`).
- 트리거: `resources_default.select=1:ncpus=32:mem=120gb`, Torque `nodes=1:ppn=16`, GPU resource custom naming.
- 제안: `select=`/`place=` parser 추가, Torque용 `nodes=...:ppn=...` 전용 parser 추가.
- 심각도: **P0 silently-wrong**.

### 5. LSF `PROCLIMIT`를 `default_cores`로 쓰는 건 의미가 틀림
- 관련 fix: **#5 미포함; 기존 parser 설계 문제** (`parsers.py:221-226`).
- 트리거: 거의 모든 LSF 9/10 production queue.
- 제안: `bqueues -l`에서 core를 직접 추론하지 말고 `bhosts -l`, host model, 또는 `LSB_DEFAULTPROCESSOR` 류 별도 probe로 분리.
- 심각도: **P0 silently-wrong**.

### 6. LSF walltime parser가 `24.0 h` 외 변형에 약함
- 관련 fix: **#5 미포함** (`parsers.py:411-418`).
- 트리거: `RUNLIMIT = 1440 min`, `7200 sec`, locale/wording 변형.
- 제안: 분/초/`HH:MM`/`HH:MM:SS`를 모두 해석하고, 가능하면 `bqueues -o ... -json` parser를 우선 사용.
- 심각도: **P0 silently-wrong**.

### 7. module false negative가 “software 없음”으로 조용히 전파됨
- 관련 fix: **#5가 scheduler 쪽만 자동화했고 software probe는 그대로** (`software.py:42-67`).
- 트리거: env-modules/Lmod init이 `.bashrc` 또는 tcsh 전용인 사이트.
- 제안: `modulecmd`/init script 존재를 직접 probe하고, bash init script를 서브셸에서 안전하게 source하는 경로를 별도 설계.
- 심각도: **P0 silently-wrong**.

### 8. Spack-only 사이트는 완전히 blind spot
- 관련 fix: **#5 미포함**; probe catalog에 Spack 관련 spec가 없다 (`probe.py:50-208`).
- 트리거: Spack shell setup + `spack load gaussian/orca` 운영 사이트.
- 제안: `spack --version`, `spack find --format`, `spack location -i` probe 추가; module보다 낮은 우선순위 fallback으로 사용.
- 심각도: **P0 silently-wrong**.

### 9. remote mode B에는 `scheduler_env` 보강이 아예 적용되지 않음
- 관련 fix: **#5의 실제 범위 한계**. local은 `env=build_scheduler_env()`지만 ssh는 raw `bash -lc`뿐이다 (`probe.py:227-234`, `248-267`).
- 트리거: workstation → login node SSH, remote node에서 scheduler path가 로그인 프로파일에만 있거나 비표준 경로인 경우.
- 제안: remote side에도 scheduler root/bin candidate bootstrap 스크립트를 넣거나, 최소한 remote `printenv` + candidate path probing을 별도 구현.
- 심각도: **P1 hard-fail**.

### 10. SGE/UGE path 후보군이 Bright/xCAT 현실을 덜 반영
- 관련 fix: **#5 신규 추가분의 path list 부족** (`scheduler_env.py:28-39`).
- 트리거: `/cm/local/apps/uge/current`, `/cm/local/apps/sge/current`, `/cm/shared/apps/uge/current`, `/cm/shared/apps/sge/current`, `/opt/gridengine`.
- 제안: trivial path-list expansion부터 하고, 장기적으로는 `SGE_ROOT`를 env + filesystem search + `qconf -help` 조합으로 찾기.
- 심각도: **P1 hard-fail**.

### 11. SLURM path 후보군도 Bright/OHPC/current symlink를 충분히 못 봄
- 관련 fix: **#5 신규 추가분의 path list 부족** (`scheduler_env.py:41-48`).
- 트리거: `/cm/local/apps/slurm/current/bin`, `/opt/slurm/current/bin`, `/opt/ohpc/pub/apps/slurm/current/bin`.
- 제안: 후보군 확장 + `scontrol`/`sinfo`가 각각 다른 디렉터리에 있는 설치를 위해 `PATH` dedupe 로직 유지.
- 심각도: **P1 hard-fail**.

### 12. PBS path 후보군이 `/opt/pbs/default/bin`과 Bright current path를 놓침
- 관련 fix: **#5 신규 추가분의 path list 부족** (`scheduler_env.py:50-57`).
- 트리거: OpenPBS/PBS Pro default symlink, `/cm/local/apps/pbspro/current/bin`, `/var/spool/pbs/bin`.
- 제안: path list 확장과 함께 `PBS_EXEC` env probe를 먼저 보도록 보강.
- 심각도: **P1 hard-fail**.

### 13. LSF는 auto-discovery/env bootstrap이 전무함
- 관련 fix: **#5 범위 밖**. `scheduler_env.py`가 LSF를 전혀 모른다.
- 트리거: `/opt/ibm/lsfsuite/lsf/<ver>/<arch>/bin`이 기본 PATH에 없는 표준 IBM 설치.
- 제안: `lsid`/`bqueues` probe, `LSF_ENVDIR`/`LSF_SERVERDIR`/`LSF_BINDIR` env probe, 대표 vendor path search 추가.
- 심각도: **P1 hard-fail**.

### 14. Cray PALS/ALPS 전용 probe가 없다
- 관련 fix: **#5 범위 밖**.
- 트리거: Cray EX/XC login node, PALS/ALPS wrapper 환경.
- 제안: PALS는 underlying SLURM와 분리해 `/opt/cray/pe/pals/default/bin` 존재 + launcher facts를 수집하고, ALPS는 사실상 unsupported로 명시하거나 별도 survey family로 분기.
- 심각도: **P1 hard-fail**.

### 15. `bash -lc` 고정은 tcsh/zsh/.bashrc-only init 사이트에서 깨진다
- 관련 fix: **#5, #6, #7 모두 bash login shell을 전제** (`probe.py:227-259`).
- 트리거: SLES tcsh default account, user가 module init을 `.cshrc`/`.zprofile`에만 둔 경우, restricted shell.
- 제안: shell detection 후 `bash`/`zsh`/`tcsh` 분기 또는 최소한 “bash login init only”를 명시하고 verify 단계에서 경고.
- 심각도: **P1 hard-fail**.

### 16. LSF JSON probe는 쏘면서도 결과를 안 쓴다
- 관련 fix: **#5가 catalog만 열어 두고 parser는 안 붙임** (`probe.py:108-121`, `survey.py:137-156`).
- 트리거: `bqueues -l` 레이아웃이 site-localized/variant라 text parse가 실패하지만 `-json`은 되는 사이트.
- 제안: `bqueues -o ... -json` parser를 먼저 붙이고 `-l`을 마지막 fallback으로 내려라.
- 심각도: **P1 hard-fail**.

### 17. mixed-scheduler gateway node에서 지금은 바로 abort
- 관련 fix: **#5가 family probe를 모두 병렬적 의미로 추가했지만 arbitration이 없다** (`survey.py:44-68`).
- 트리거: migration 중인 site, admin/gateway node에 `sinfo`와 `qstat`가 함께 있는 경우.
- 제안: topology evidence, env namespace (`SLURM_`, `PBS_`, `SGE_`, `LSB_`), active daemon reachability를 가중치로 ranking.
- 심각도: **P1 hard-fail**.

### 18. `build_scheduler_env()`가 1회 캐시라 late export/module load를 못 본다
- 관련 fix: **#5가 LRU cache를 도입** (`scheduler_env.py:243-255`).
- 트리거: 첫 probe 실패 후 같은 long-lived process에서 `export SGE_ROOT=...` 또는 `module load slurm`를 한 경우.
- 제안: refresh/verify 경로에서 cache bust hook 추가, 또는 env fingerprint 기반 cache key로 변경.
- 심각도: **P1 hard-fail**.

### 19. `qhost` fallback이 cluster-global min이라 queue별로 너무 보수적일 수 있음
- 관련 fix: **#6/#7 도입 로직의 부작용** (`parsers.py:235-263`, `survey.py:173-217`).
- 트리거: fat queue + serial queue + GPU queue가 한 클러스터에 공존하는 mixed-node site.
- 제안: `qhost -q` 또는 queue-instance mapping으로 queue별 min resource를 계산.
- 심각도: **P2 cosmetic / conservative-underestimate**.

### 20. `_find_pbs_bin()`의 `qstat -Q` 휴리스틱은 다소 느슨함
- 관련 fix: **#5 신규 path bootstrap의 판별 로직** (`scheduler_env.py:207-224`).
- 트리거: 비표준 `qstat` wrapper가 `-Q`에 0을 반환하는 특수 site.
- 제안: stdout pattern을 더 엄격히 확인하고, `qstat --version` 또는 `pbsnodes` 존재를 교차 검증.
- 심각도: **P2 cosmetic / rare-misclassification**.

### 21. `module -t avail` normalization은 노이즈를 충분히 못 걷어냄
- 관련 fix: **#5 미포함; software survey 기존 한계** (`software.py:215-233`).
- 트리거: `(D)` marker, hidden module note, localized header, `.lua` suffix noise, category prefix.
- 제안: Lmod/env-modules 별 normalization fixture를 추가하고, hidden/default marker stripping을 강화.
- 심각도: **P2 cosmetic**.

정리:

- **P0 = 8건**: 1-8
- **P1 = 10건**: 9-18
- **P2 = 3건**: 19-21

## 6. Recommended next waves (prioritized backlog)

### (a) trivial path-list additions closing many holes cheap

1. `scheduler_env.py` 후보군 즉시 확장.
   - SGE/UGE: `/cm/local/apps/uge/current`, `/cm/local/apps/sge/current`, `/cm/shared/apps/uge/current`, `/cm/shared/apps/sge/current`, `/opt/gridengine`.
   - SLURM: `/cm/local/apps/slurm/current/bin`, `/opt/slurm/current/bin`, `/opt/ohpc/pub/apps/slurm/current/bin`.
   - PBS: `/opt/pbs/default/bin`, `/var/spool/pbs/bin`, `/cm/local/apps/pbspro/current/bin`, `/cm/shared/apps/pbspro/current/bin`.
2. LSF는 path list만으로도 초반 효과가 크다.
   - `/opt/ibm/lsfsuite/lsf/*/bin`, `/cm/shared/apps/lsf/current/*/bin` 정도만 추가해도 체감 차이가 큼.
3. docs wording도 함께 좁혀야 한다.
   - `docs/installation-troubleshooting.md:246-252`의 문구는 “표준 SGE/SLURM/PBS local install에서는 보통”으로 범위를 줄여야 안전하다.

### (b) parser hardening for variant outputs

1. **LSF JSON parser를 최우선**으로 붙여라.
   - 지금 가장 값싸게 많은 구멍을 막는 작업이다.
2. **PBS `select=` / Torque `nodes:ppn` parser**를 분리해라.
   - OpenPBS/PBS Pro/Torque를 한 parser로 뭉개는 건 장기적으로 유지보수 비용이 더 크다.
3. **SGE state/slots semantics**를 손봐라.
   - `slots<=1` fallback만으로는 production UGE/SoGE를 커버하지 못한다.
4. **SLURM version-skew fixture**를 늘려라.
   - 최소 `sinfo --json` 실패 + `scontrol --oneliner`만 성공하는 18/20 시나리오와, `DefMemPerCPU`/`TRES` 중심 22/23 시나리오를 분리.

### (c) module/shell-init compat

1. `which module` 의존을 낮추고 `modulecmd`/init script existence probe를 추가.
2. shell 종류를 먼저 식별하고, 적어도 `bash`와 `tcsh`는 분기.
3. `ml` alias를 직접 probe할 필요는 낮지만, **“module function 없음 = module stack 없음”**이라는 결론은 금지.
4. Spack probe를 추가해 module-less cluster도 software survey에 걸리게 하라.
5. remote mode B에도 local mode A와 동등한 bootstrap 전략을 설계하라.

### (d) deeper architecture

1. **scheduler family별 topology probe**를 더 깊게 나눠라.
   - LSF: `lsid`, `bqueues`, 필요 시 `bhosts`.
   - PBS: `qstat`, `pbsnodes`, `qmgr -c 'p s'`까지는 아니더라도 queue/exec facts 분리.
   - SGE/UGE: `qhost -q`, `qconf -se`, hostgroup mapping.
2. **queue facts와 host facts를 분리**해라.
   - 지금은 queue와 노드 정보를 한 구조체에 억지로 합쳐서 mixed-node cluster에서 오류가 난다.
3. **remote bootstrap parity**를 만들라.
   - local만 `build_scheduler_env()`를 쓰는 현재 구조는 기능적으로 반쪽이다.
4. **arbitration layer**를 넣어 mixed-scheduler node에서도 best-effort 선택이 되게 하라.

## 7. Out-of-scope but flag

- **container-only login nodes / Apptainer(Singularity) wrapper**
  - host의 scheduler/module init이 container 내부 `bash -lc`에 어떻게 보이는지 현재 모델은 전혀 설명하지 못한다.
- **Bright cluster manager paths**
  - 이번 감사에서 가장 반복적으로 나온 path hole이 바로 `/cm/local/apps/*/current`였다.
- **OOD / JupyterHub / VS Code Server launch context**
  - non-login shell, sanitized env, web session proxy 때문에 local mode A와 실제 사용자 체감이 달라질 가능성이 크다.
- **xCAT / Cray PE shell bootstrap**
  - site wrapper가 `/etc/profile.d` 외부에 있거나 tcsh 전용이면 지금 probe 구조와 충돌한다.

## Bottom line

PR #81은 **“SGE 8.1.8 / CentOS 7 / bash login shell / local mode A”에서는 의미 있는 전진**이다. 다만 그 성공을 **LSF, Torque, Bright-managed UGE, Cray, Spack-only, tcsh-default** 환경으로 일반화하면 안 된다. 다음 웨이브는 새 기능보다도 **경로 후보군 확장 + parser variant hardening + module/shell bootstrap 재설계**가 먼저다.
