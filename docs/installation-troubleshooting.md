# chemsmart 설치 오류 보고서

**작성일:** 2026-05-11  
**환경:** CentOS 7 (Linux 3.10.0), Python 3.10, conda (miniconda3), SGE 8.1.8 HPC 클러스터 헤드 노드

---

## 요약

`environment.yml` 기반으로 `chemsmart` conda 환경을 신규 구성하는 과정에서 총 **6가지** 오류가 발견되었으며, 모두 수정 완료하였다.

---

## 오류 1: `rdkit==2025.3.3` — pip에서 Python 3.10 미지원

### 증상
```
ERROR: Could not find a version that satisfies the requirement rdkit==2025.3.3
ERROR: No matching distribution found for rdkit==2025.3.3
```

### 원인
`pyproject.toml`이 `rdkit==2025.3.3`으로 핀고정되어 있으나, PyPI의 rdkit 휠은 Python 3.10 환경에서 최대 `2024.3.2`까지만 배포된다. rdkit 2025.x 이상의 휠은 Python 3.11+ 전용이다.

### 해결
- conda-forge에서 최신 rdkit(`2026.03.2`)을 설치한 뒤, pyproject.toml의 버전 핀을 `rdkit>=2024`로 완화하여 pip가 conda 설치본을 유효한 것으로 인식하도록 수정.

```diff
# pyproject.toml
-    "rdkit==2025.3.3",
+    "rdkit>=2024",
```

```bash
conda install -n chemsmart -c conda-forge rdkit
```

---

## 오류 2: `spglib`, `pyzmq` 소스 빌드 실패 — 구형 시스템 GCC

### 증상
```
error: 'for' loop initial declarations are only allowed in C99 mode
ninja: build stopped: subcommand failed.
× Failed to build installable wheels for some pyproject.toml based projects
╰─> spglib, pyzmq
```

### 원인
서버 OS가 CentOS 7이며, 시스템 GCC(4.8 계열)의 기본 C 표준이 C89이다. `spglib`와 `pyzmq`는 C99 이상을 요구하는 코드를 포함하고 있어 빌드 실패한다. PyPI에 Python 3.10용 manylinux 바이너리 휠이 없는 버전이 요청된 경우 소스 빌드를 시도한다.

### 해결
두 패키지를 conda-forge에서 선빌드된 바이너리로 설치하면 컴파일 없이 해결된다.

```bash
conda install -n chemsmart -c conda-forge spglib pyzmq
```

---

## 오류 3: `datetime.UTC` — Python 3.10 미지원 심볼

### 증상
```
ImportError: cannot import name 'UTC' from 'datetime'
```

`chemsmart agent` 서브커맨드 실행 시 발생하여, agent 그룹 전체가 `_missing_agent_group`으로 대체되어 `No such command 'doctor'` 오류로 이어진다.

### 원인
`datetime.UTC` 상수는 Python 3.11에서 추가되었다. `environment.yml`은 Python 3.10을 지정하고 있어 해당 심볼을 임포트할 수 없다.

### 영향 파일 (소스 5개 + 테스트 3개)

| 파일 | 위치 |
|------|------|
| `chemsmart/agent/cli.py` | 소스 |
| `chemsmart/agent/core.py` | 소스 |
| `chemsmart/agent/tui/services/job_poller.py` | 소스 |
| `chemsmart/agent/wizard/refresh.py` | 소스 |
| `chemsmart/agent/wizard/cache.py` | 소스 |
| `tests/agent/wizard/test_cache.py` | 테스트 |
| `tests/agent/wizard/test_refresh.py` | 테스트 |
| `tests/agent/test_agent_cli_run.py` | 테스트 |

### 해결

모든 파일에서 아래와 같이 수정:

```diff
-from datetime import UTC, datetime
+from datetime import datetime, timezone
+
+UTC = timezone.utc
```

(`cache.py`처럼 `timedelta`도 함께 import하는 경우는 기존 import line에 통합)

---

## 오류 4: `PydanticSchemaGenerationError` 미처리 — `ToolRegistry` 빌드 실패

### 증상
```
pydantic.errors.PydanticSchemaGenerationError: Unable to generate pydantic-core schema
for <class 'chemsmart.io.molecules.structure.Molecule'>
```

`chemsmart agent doctor` 실행 시 traceback과 함께 종료된다.

### 원인
`chemsmart/agent/registry.py`의 `_schema_friendly_annotation` 함수가 Pydantic v2의 `PydanticSchemaGenerationError`를 catch하지 않는다. `Molecule`, `Job`, `SubmitTransport` 등 Pydantic 모델이 아닌 도메인 클래스가 tool 시그니처에 등장할 때 이 예외가 발생한다. 오류가 발생하면 `ToolRegistry.default()` 전체가 실패하여 14개 tool이 등록되지 않는다.

### 해결

```diff
# chemsmart/agent/registry.py
-from pydantic.errors import PydanticInvalidForJsonSchema
+from pydantic.errors import PydanticInvalidForJsonSchema, PydanticSchemaGenerationError

 def _schema_friendly_annotation(...):
     try:
         TypeAdapter(annotation).json_schema()
-    except (PydanticInvalidForJsonSchema, TypeError):
+    except (PydanticInvalidForJsonSchema, PydanticSchemaGenerationError, TypeError):
```

해당 타입은 `Any`로 폴백 처리되며 tool 동작에는 영향 없다.

---

## 오류 5: SGE 큐 이름의 점(`.`) — `validate_identifier` 거부

### 증상
```
chemsmart.agent.wizard.probe.ProbeError: Invalid identifier slot value: '20core.q'
```

`chemsmart agent wizard <name>` 실행 시, `qconf -sql`이 반환한 큐 이름(`20core.q`, `40core.q`, `fill_up`)을 처리하는 과정에서 발생한다.

### 원인
`chemsmart/agent/wizard/probe.py`의 identifier 검증 패턴이 점(`.`)을 허용하지 않는다.

```python
_IDENTIFIER_PATTERN = re.compile(r"^[\w\-]+$")  # 점 불허
```

SGE/GridEngine의 큐 이름 관례는 `<name>.q` 형식으로 점을 포함한다.

### 해결

```diff
-_IDENTIFIER_PATTERN = re.compile(r"^[\w\-]+$")
+_IDENTIFIER_PATTERN = re.compile(r"^[\w\-\.]+$")
```

---

## 오류 6: SGE `h_vmem=INFINITY` — `MEM_GB` null로 검증 실패

### 증상
```
Error: wizard output did not validate:
- Missing required SERVER.MEM_GB.
```

`wizard` 명령이 YAML을 생성하지만, `MEM_GB: null`이어서 자체 검증에서 실패한다.

### 원인
이 서버의 SGE 큐는 `h_vmem=INFINITY`로 설정되어 있어(메모리 상한 없음), `parse_sge_qconf_sq`의 메모리 파서가 `None`을 반환한다. `wizard`는 `qhost` 명령으로 실제 노드 메모리를 조회하는 fallback 로직이 없었다.

실제 노드 구성:
- node02–node06: 20코어, 메모리 62.7 GB
- node07–: 40코어, 메모리 187.4 GB

### 해결

3개 파일을 수정하여 `qhost` 출력을 fallback 메모리 소스로 사용하도록 구현:

**`probe.py`** — `survey.sge.qhost` probe spec 추가:
```python
"survey.sge.qhost": ProbeSpec(
    template_id="survey.sge.qhost",
    argv_template=("qhost",),
    slot_validators={},
),
```

**`parsers.py`** — `parse_sge_qhost()` 함수 추가 (각 노드 `MEMTOT` 컬럼 파싱, 최솟값 반환):
```python
def parse_sge_qhost(payload: str) -> int | None:
    mem_values: list[int] = []
    for line in payload.splitlines():
        parts = line.split()
        if len(parts) < 8:
            continue
        hostname = parts[0]
        if hostname in ("HOSTNAME", "global", "-"):
            continue
        parsed = _parse_mem_gb(parts[7])
        if parsed is not None:
            mem_values.append(parsed)
    return min(mem_values) if mem_values else None
```

**`survey.py`** — `_probe_sge()`에 fallback 로직 추가:
- `qhost` 1회 호출 → 노드 최소 메모리 취득
- 각 큐의 `default_mem_gb is None`인 경우 해당 값으로 보완

### 결과

`MEM_GB: 63`(= 62.7 GB 올림)으로 YAML 생성 및 검증 통과:

```yaml
SERVER:
  SCHEDULER: SGE
  QUEUE_NAME: all.q
  NUM_HOURS: 24
  MEM_GB: 63
  NUM_CORES: 1
  ...
ORCA:
  EXEFOLDER: /usr/bin
  ...
```

---

## 최종 테스트 결과

| 명령어 | 결과 |
|--------|------|
| `chemsmart agent doctor` | ✅ ping ok (gpt-5.4-2026-03-05) |
| `chemsmart agent tools` | ✅ 14개 tool 등록 |
| `chemsmart agent sessions` | ✅ 정상 |
| `chemsmart agent ask "..."` | ✅ LLM 응답 스트리밍 정상 |
| `chemsmart agent run ... --dry-submit` | ✅ Gaussian 입력 파일 생성 |
| `chemsmart agent wizard mycluster` | ✅ SGE 감지 + YAML 생성 (MEM_GB: 63) |
| `chemsmart agent wizard-verify <name>` | ✅ 정상 |
| `chemsmart agent wizard-refresh <name>` | ✅ 정상 |
| wizard 유닛 테스트 (78개) | ✅ 전부 통과 |

---

## 참고: `wizard` 실행 전 수동 scheduler 로드는 보통 필요 없음

이제 `scheduler_env.py`가 표준 설치 경로에서 SGE/SLURM/PBS 환경을 자동 탐지한다.
`/opt/sge`, `/opt/uge`, `/cm/shared/apps/sge` 같은 일반적인 위치라면
`source .../settings.sh`를 먼저 실행하지 않아도 `chemsmart agent wizard`가 바로 동작한다.
다만 SGE가 비표준 경로에 설치된 환경이라면, wizard 실행 전에
`export SGE_ROOT=/custom/path/to/sge`처럼 `SGE_ROOT`를 먼저 지정하면 된다.
