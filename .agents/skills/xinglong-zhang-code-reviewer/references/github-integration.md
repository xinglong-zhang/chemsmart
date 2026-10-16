# GitHub integration

This skill supplies review guidance; it does not configure GitHub or enable automatic reviews by itself. Install the complete skill directory, including its `references/` and `assets/` resources, in a location available to the agent that performs the review.

## Repository setup

1. Add the skill to the repository at `.agents/skills/xinglong-zhang-code-reviewer/`, preserving its directory structure.
2. Merge the relevant section from [github-AGENTS.md](../assets/github-AGENTS.md) into the repository's existing `AGENTS.md` or agent instructions. Keep unrelated instructions and do not replace an existing file.
3. Confirm that the agent can read the skill and referenced resources in the environment where reviews run.

## Manual reviews

Ask the configured agent to review a specific pull request and provide its URL or number and any scope constraints. The agent should inspect the PR diff and relevant code, then return actionable findings with evidence and validation limitations. A request to review does not by itself authorize posting a review, changing the branch, or merging the PR.

## Automatic reviews

Automatic review requires a GitHub-supported review integration and its required repository or organization settings and permissions. Configure that integration separately, following the current GitHub documentation for the chosen feature. Make the skill available to the configured agent and add repository instructions that direct it to use the skill. Verify the integration with a test PR; installing the skill alone does not create a workflow or trigger.
