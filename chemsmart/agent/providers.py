"""One declarative registry of agent providers.

Adding a provider means: one declaration here, its adapter module under
``chemsmart/agent/runtime/``, and the two protocol construction sites
(``runtime_config`` and the unified session). Every other surface -- the
``chemsmart config agent`` menus, the default key labels and engine-env
scrub tokens, the profile loader's endpoint recognition and effort
vocabulary, and the loop's runnable set -- derives from this table, so
"this provider exists / is runnable" is stated exactly once.
"""

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping


@dataclass(frozen=True)
class ProviderDeclarationV1:
    """Everything the non-protocol surfaces need to know about a provider."""

    name: str
    display_name: str
    endpoint: str
    key_labels: tuple[str, ...]
    key_label_token: str
    #: Validated effort vocabulary; an empty tuple accepts any value.
    efforts: tuple[str, ...]
    default_effort: str
    effort_error: str
    #: What ``chemsmart config agent`` offers ("" renders as "(omit)").
    cli_efforts: tuple[str, ...]
    #: Help text only -- the model prompt has NO preselection.
    model_hint: str
    requires_preserved_thinking: bool = False
    runnable: bool = True
    refusal: str = ""


PROVIDERS: Mapping[str, ProviderDeclarationV1] = MappingProxyType(
    {
        "openai": ProviderDeclarationV1(
            name="openai",
            display_name="OpenAI",
            endpoint="https://api.openai.com/v1",
            key_labels=("OPENAI_API_KEY",),
            key_label_token="OPENAI",
            efforts=("", "low", "medium", "high"),
            default_effort="",
            effort_error=(
                "OpenAI reasoning effort must be low, medium, high, or "
                "omitted"
            ),
            cli_efforts=("", "low", "medium", "high"),
            model_hint="the exact model id from your OpenAI account",
        ),
        "anthropic": ProviderDeclarationV1(
            name="anthropic",
            display_name="Anthropic",
            endpoint="https://api.anthropic.com",
            key_labels=("ANTHROPIC_API_KEY",),
            key_label_token="ANTHROPIC",
            efforts=(),
            default_effort="",
            effort_error="",
            cli_efforts=("", "low", "medium", "high"),
            model_hint="the exact model id from your Anthropic account",
            runnable=False,
            refusal=(
                "the anthropic adapter is not registered in this release; "
                "the profile remains valid configuration"
            ),
        ),
        "alibaba-token-plan": ProviderDeclarationV1(
            name="alibaba-token-plan",
            display_name="Alibaba Token Plan",
            endpoint=(
                "https://token-plan.ap-southeast-1.maas.aliyuncs.com"
                "/compatible-mode/v1"
            ),
            key_labels=("ALIBABA_TOKEN_PLAN_KEY",),
            key_label_token="ALIBABA",
            efforts=("high", "max", "xhigh"),
            default_effort="max",
            effort_error=(
                "Alibaba Token Plan reasoning effort must be high, max or "
                "xhigh"
            ),
            cli_efforts=("high", "max", "xhigh"),
            model_hint="the exact catalog id your Token Plan serves",
            requires_preserved_thinking=True,
        ),
        "deepseek": ProviderDeclarationV1(
            name="deepseek",
            display_name="DeepSeek",
            endpoint="https://api.deepseek.com",
            key_labels=(
                "DEEPSEEK-api-key",
                "DEEPSEEK_API_KEY",
                "DEEPSEEK-API-KEY",
            ),
            key_label_token="DEEPSEEK",
            efforts=("high", "max"),
            default_effort="max",
            effort_error="DeepSeek reasoning effort must be high or max",
            cli_efforts=("high", "max"),
            model_hint="the exact model id from your DeepSeek account",
        ),
    }
)


def declaration_for_endpoint(endpoint: str) -> ProviderDeclarationV1 | None:
    for declaration in PROVIDERS.values():
        if declaration.endpoint == endpoint:
            return declaration
    return None


def runnable_provider_names() -> frozenset[str]:
    return frozenset(
        name for name, declaration in PROVIDERS.items() if declaration.runnable
    )
