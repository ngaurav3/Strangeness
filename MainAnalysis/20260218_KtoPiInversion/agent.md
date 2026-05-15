# Agent Digital Twin

This file is a practical digital-twin profile of the user for Codex in this
repository. It is not a generic style guide. It reflects how the user actually
drives DELPHI analysis work in this project.

## 1. Core Working Pattern

The user solves problems iteratively and operationally.

Typical pattern:

1. pick a concrete target
2. make a direct change
3. rebuild or rerun
4. inspect the output
5. adjust details quickly
6. push when satisfied

This means Codex should optimize for:

- momentum
- correctness
- explicit verification
- low-friction iteration

Do not default to long design discussion when the user is clearly in execution
mode.

## 2. What the User Cares About

The user consistently values:

- code-grounded answers
- exact observable definitions
- figure correctness
- note/paper synchronization
- reproducible documentation
- method validation before claiming a workflow change

The user is comfortable changing strategy midstream if the validation argues for
it.

## 3. How the User Breaks Down Problems

The user tends to decompose work into three layers:

### Physics layer

- What is the observable?
- What is included or excluded?
- Is the interpretation defensible?

### Analysis layer

- What is the active workflow?
- Is the closure/refolding/systematics behavior acceptable?
- Is the chosen method actually better, or just visually nicer?

### Presentation layer

- Does the figure say the right thing?
- Is the legend clean?
- Does the note wording match the code?
- Is the paper consistent with the note?

Codex should usually follow the same decomposition.

## 4. Default Response Strategy for This User

When the user asks a technical question:

- audit the active implementation first
- answer directly
- then state the caveat that matters

When the user asks for a change:

- make the change
- rebuild/regenerate the affected artifact
- verify the result
- push only when asked, or when the request clearly implies it

When the user asks for a methodological comparison:

- compare closure
- compare refolding
- compare systematic impact
- compare interpretability
- make a decision, not just a list

## 5. User Preferences in Communication

The user prefers:

- short factual updates while work is in progress
- final answers that start with the conclusion
- explicit file references for important changes
- no fluff
- no unnecessary motivational language

The user does not want vague statements like:

- "it should be fine"
- "probably"
- "I think this is okay"

If there is uncertainty, quantify it or point to the exact unresolved piece.

## 6. Validation Standard

The user repeatedly pushes on one question:

- does the analysis actually validate the claimed method?

For method changes, Codex should not stop at:

- "the new plot looks better"

Codex should check:

- closure
- refolding
- overflow/tail behavior
- systematic comparison
- whether the workflow is physically consistent

If a new method is better only in some bins but worse globally, say that
clearly.

## 7. Repository Behavior Expected by the User

The user expects Codex to understand the repo split:

1. main code repo
2. analysis-note repo
3. paper-draft repo
4. external exploratory repos

Operational rule:

- do not mix these scopes casually
- keep commits targeted
- know when a push request applies to the note only, the paper only, or the code
  repo

## 8. Scientific Caution Areas the User Repeatedly Returns To

These topics are persistent high-sensitivity areas:

- exact event selection
- fiducial versus whole-event interpretation
- weak-decay-daughter veto
- heavy-flavor feed-down
- PID scale-factor convention
- generator settings, tunes, and versions
- current workflow versus legacy workflow
- whether a cross-check is strong enough to replace the nominal method

Codex should surface these caveats proactively when relevant.

## 9. Figure and Note Workflow Expectations

The user often works through figures as the practical interface to the physics.

That means:

- if a figure changes, make sure the underlying method is still correct
- if the method changes, make sure the figure and caption are updated
- if the note changes, check whether the paper should also be updated

The user also often refers to figures by number. Codex should verify the actual
current figure identity when numbering may have shifted.

## 10. Documentation Expectation

The user wants a written trace of significant work.

Default behavior:

- add a markdown note in `report/` for nontrivial studies
- update `worklog.md` for user-prompt history
- keep the main technical summary current in `Delphi.md`

If a conclusion required real investigation, leave a written record.

## 11. Decision Heuristics for Codex

When deciding how to act, use these defaults:

- prefer the active validated workflow over a newer but weaker one
- prefer code-backed answers over note-backed answers if they disagree
- prefer a documented caveat over silent ambiguity
- prefer a targeted commit over a broad commit
- prefer an explicit cross-check over hand-waving

## 12. Short Operational Summary

If Codex has to compress the user model to one paragraph, use this:

The user is an iterative physics-analysis operator who wants rapid execution,
hard validation, exact observable definitions, and synchronized documentation.
They are willing to move fast, but only if closure, refolding, and the code path
support the claim. Codex should act directly, verify concretely, document
nontrivial conclusions, and keep the main code repo, note repo, and paper repo
cleanly separated.
