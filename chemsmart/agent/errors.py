"""Public error types of the agent layer.

`ContractError` is the one exception every agent surface raises for a
violated contract; outside layers import it from here rather than from
the private ``_contracts`` module.
"""

from chemsmart.agent._contracts import ContractError

__all__ = ["ContractError"]
