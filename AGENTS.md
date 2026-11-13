# CHEMSMART agent instructions

## Code reviews

For code reviews, read and apply:
`.agents/skills/xinglong-zhang-code-reviewer/SKILL.md`

Separate correctness defects, OOP/API recommendations, and style preferences.
Support findings with inspected code and concrete consequences.

If the skill files are unavailable, report that limitation rather than
claiming to have applied them.

## Improving the code-review skill

After Xinglong corrects a review, propose the corresponding skill improvement.

When Xinglong explicitly says “learn this”:
- Read the current skill before editing.
- Update the relevant rule in SKILL.md.
- Record the feedback and its source in references/review-evidence.md.
- Add a regression example in references/review-cases.md, creating it
  if necessary.
- Validate the change and prepare a PR explaining the improvement.

Give explicit corrections precedence over older guidance.
Do not modify the skill merely because a PR author disagrees with a finding.
