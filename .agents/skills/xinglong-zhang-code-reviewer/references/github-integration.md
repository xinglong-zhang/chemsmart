# GitHub review integration

Use this reference when installing the reviewer in a repository or
explaining PR triggers. Recheck official documentation before changing
product-specific setup.

## Store the review standards with the code

Store the complete skill directory at:

.agents/skills/xinglong-zhang-code-reviewer/

Include SKILL.md, agents/openai.yaml, and all referenced resources.
Keep relative paths intact.

Merge assets/github-AGENTS.md into the repository root AGENTS.md when
setting up review instructions. Preserve existing instructions,
including the feedback and skill-improvement workflow.

A personal ChatGPT skill is not automatically available to a GitHub
reviewer. The repository copy provides the shared review standards.

## Enable GitHub reviews

Connect the repository to Codex and enable code review for that
repository in Codex settings.

To request a review, add an ordinary PR conversation comment:

```text
@codex review
```

To provide a focus:

```text
@codex review for scientific correctness and object-ownership/API risks,
following AGENTS.md.
```

For automatic reviews of students' PRs, enable repository automatic
review for all PRs and check the configured triggers and personal
preferences.

Adding skill files does not itself enable the GitHub integration.

## Request a full assessment

The GitHub review host may restrict which priorities it publishes.
Respect its current reporting threshold. Never increase a finding's
severity merely to make it appear.

Use a normal Codex task or local session for a full assessment covering
correctness, OOP/API recommendations, and style preferences.

For local Codex:

```text
Use $xinglong-zhang-code-reviewer to review this PR against its base
branch. Include correctness findings, OOP/API recommendations, and
separate style notes.
```

For a separate GitHub Codex task, explicitly identify the instruction
file:

```text
@codex Read .agents/skills/xinglong-zhang-code-reviewer/SKILL.md and
assess this PR. Include correctness findings, OOP/API recommendations,
and separate style notes. Return the assessment without changing code.
```

Treat this as a task prompt, not an independently tested integration.
Do not assume that $skill-name alone is a GitHub comment command.

## Verify the installation

- Confirm that all relative resources referenced by SKILL.md exist.
- Confirm that icon paths in agents/openai.yaml resolve from the
  skill directory.
- Confirm that root AGENTS.md points to the correct skill path.
- Check that reviews cite the inspected PR revision.
- Distinguish file validation from a successful live GitHub review.

For review-quality checks, examine an object-owned operation, a
legitimate pure helper, and an unrelated change. Confirm useful
ownership feedback without forcing every function into a class.

## Troubleshooting

If a PR comment receives no response:

- Confirm that the PR targets the connected repository.
- Confirm that code review is enabled for that repository.
- Test the exact trigger @codex review.
- For automatic reviews, inspect repository settings, personal
  preferences, and the configured PR event triggers.

If a resource is missing, restore it at the referenced relative path.
Do not claim to have read unavailable instructions.

## Official sources

- GitHub integration:
  https://learn.chatgpt.com/docs/third-party/github
- Repository skills:
  https://learn.chatgpt.com/docs/build-skills
- Custom review rules:
  https://developers.openai.com/blog/custom-code-review-rules-for-codex
