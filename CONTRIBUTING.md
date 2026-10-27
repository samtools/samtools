Contribution guide
==================

We value all contributions that lead to a net improvement in the
quality of this project.  These include but are not limited to

- Requests for new functionality
- Implementations of new functionality (but please ask first)
- Reports of bugs
- Bug fixes
- Documentation (issues and/or improvements)
- Performance, robustness and security improvements (contact us privately)

However we must ensure for all changes that copyright is compatible
with our license (see LICENSE file).  We also have requirements of git
commit styles and AI policies.  If you are contributing changes
please read this policy or your PR may be rejected.


Complexity
==========

We accept that some pieces of work are simply big by their nature, but
please do keep every pull request (PR) as small as is possible in
order to speed up review.

- Do not mix multiple pieces of functionality in a single PR unless they
  are explicitly connected and all are required in order to pass the
  tests.  Functionality changes, documentation and test cases are fine
  together within a single PR.

- Consider breaking down large PRs if possible, although do not do
  this arbitrarily if it means sub-PRs will fail tests until they are
  all merged in.

- If you have a large piece of work, create an issue before hand so we
  can discuss the best approach and scope of the changes.  This will
  make it easier for us to review and may save you wasted effort.

- PRs should be against a reasonably modern develop branch.  We will
  do some rebasing and maybe squashing during merge (see `git rebase
  -i`), but dealing with conflicts will slow down our review or may
  even cause it to be rejected.

- Break the changes down into distinct commits where the logical
  functionality is split into discrete components.

- Do not add new dependencies unless there is a very good reason.  We
  want to keep the package as simple to build as possible.

- Do not change the API or ABI without prior discussion.


Completeness
============

We require documentation (where relevant) and test cases for all PRs.

- Commit messages should be sufficient to understand the intent of the
  changes.  Please do not have a short one-line commit message unless
  the change itself is trivial.

- Do not include huge test files.  If it requires a lot of data for
  validation, consider auto-generating it (in a deterministic way).
  Do not attempt to access remote sources during testing.

- Ensure you only use public data or data that you have permission to
  publish (see "Signing" below).

- The CI system performs some basic code validity checks, so ensure
  these pass before making a PR.  We will also check for correct
  handling of error cases and recovery.


Signing your work
=================

Every PR must be explicitly stated by you to adhere to our policies.
We follow the Developer Certificate of Origin (DCO) policy produced by
the Linux Foundation: See https://developercertificate.org/
A copy of this is included below:

    Developer Certificate of Origin
    Version 1.1
    
    Copyright (C) 2004, 2006 The Linux Foundation and its contributors.
    
    Everyone is permitted to copy and distribute verbatim copies of this
    license document, but changing it is not allowed.
    
    
    Developer's Certificate of Origin 1.1
    
    By making a contribution to this project, I certify that:
    
    (a) The contribution was created in whole or in part by me and I
        have the right to submit it under the open source license
        indicated in the file; or
    
    (b) The contribution is based upon previous work that, to the best
        of my knowledge, is covered under an appropriate open source
        license and I have the right under that license to submit that
        work with modifications, whether created in whole or in part
        by me, under the same open source license (unless I am
        permitted to submit under a different license), as indicated
        in the file; or
    
    (c) The contribution was provided directly to me by some other
        person who certified (a), (b) or (c) and I have not modified
        it.
    
    (d) I understand and agree that this project and the contribution
        are public and that a record of the contribution (including all
        personal information I submit with it, including my sign-off) is
        maintained indefinitely and may be redistributed consistent with
        this project or the open source license(s) involved.


To certify your agreement with the above, add a Signed-off-by line to
the end of each commit message. For example:

    Signed-off-by: William Wordsworth <bill@poets.org>

Please do not use anonymous contributions or rely on GitHub username
aliases.  The sign off needs to be a real name.

You can do this automatically via `git commit -s` if you have filled
out the `[user]` section of ~/.gitconfig.

If multiple authors are involved, we need a sign off statement from
every author.


AI Policy
=========

(Credit to: https://docs.kernel.org/process/coding-assistants.html)

We acknowledge that the use of AIs and LLMs can be a productivity aid.
However it is still your responsibility to ensure all generated code
is validated for accuracy and if appropriate brevity.  We do not have
the time to review many large AI-generated PRs. (See Complexity above.)

AI agents are not permitted to add Signed-off-by tags.  We will only
accept DCO sign off by humans.

As part of the sign off, you are responsible for:

- Reviewing all code generated by an automated tool (including AIs).
  Every line needs reviewing.

- Taking full responsibility for all code and interactions during the
  review process.

- Ensuring compatibility with our software license.

- Adding the appropriate Signed-off-by tags, one per author involved.

To aid tracking of which AI tools are utilised, please add an
"Assisted-by:" line to each commit as appropriate.  These should take
the format of:

    Assisted-by: AGENT_NAME:MODEL_VERSION [TOOL ...]

Where:

- AGENT_NAME is the name of the AI (e.g. ChatGPT, Claude, Gemini)

- MODEL_VERSION is the model name and version (e.g. gpt-5.5,
  claude-3.7-sonnet, gemini-3-flash).

- TOOL is an optional list of specialised analysis tools,
  (e.g. coccinelle or smatch).

For ease of parsing, replaces spaces in any of these fields with hyphens.
We do not need information on basic tools such as compilers or non-AI
assisted editors.

For example:

    Assisted-by: Claude:claude-3.7-sonnet sparse smatch

Please also add an Assisted-by statement to the end of your pull request
description too.

Commit messages and PR descriptions should be written by a human.
