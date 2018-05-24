# Contribution guidelines

We welcome any bug-fixes and patches to this project. However, we request that
your contributions be well-tested and follow our style guidelines to maintain
the code quality of this project. Please follow these guidelines when providing
contributions.
  
## Coding guidelines

If in doubt, use [Google's
python](https://google.github.io/styleguide/pyguide.html) coding guidelines.
Some minor changes have been done.

- Ensure textwidth of 80 characters is followed.
- Always use semi-colon
- List arguments of the functions in separate lines so that the arguments
  look like a list of items.
    - If you can fit all the variables (in their respective new lines) within
      80 characters, align the variable list at the column after the `(`.
      Otherwise, break the line with `...`, add one indent, and start the
      variable list.
    - Name-value pairs appear in the same line.
    - Names are strings written in upper camel case.
    - Few functions are allowed exception to this rule for readability.
        * `validateattributes`
        * `isa`
- Naming rules:
    - Use `obj` to refer to the object within the class methods.
    - Use descriptive function names (but not too long).
        - Verbs to be used:
            - get: Will return a value described by the function
            - compute: Updates internal variables as described by the function
    - Write class names in upper camel cased.
    - Write function names in lower camel cased.
    - Write variable names in lower cased words separated by underscores.
    - Write acronyms keeps the appropriate camel case (e.g. `LTISystem` class
      becomes `LtiSystem`, and `getLTISystemMatrix()` becomes
      `getLtiSystemMatrix()`.).
- Use spaces between operators (`+,-,/,*,=`) and assignments for readability.
- Use `assert`/`error` functions to encode sanity checks into the code.
    - Leverage input handling in the existing code to avoid redundancy. However,
      mark the leveraging action by adding in the description comment `DELEGATED
      INPUT HANDLING`.
    - Avoid `if-condition-error-end` for sanity checks. Use `assert`.
    - Exceptions must include a `msgIdentifer` of the form `SReachTools:XXX`.
      Identifiers currently in use are
        * `SReachTools:internal`    --- For internal errors (TODOs).
        * `SReachTools:invalidArgs` --- For invalid arguments.
        * `SReachTools:setup_error` --- For lack of the necessary dependencies.
    - Ensure that exception message is succinct and disambiguous. Don't be TeX!
    - Test all exceptions provided in the code.
    - Exception message do not end with a period.
    - Exception messages need not be complete sentences in the interest of
      brevity.
    - Ensure that `assert`/`error` function calls have their arguments in the
      list form.
- Ensure the help tags of the classes and functions follow the same
  structure as `LtiSystem` and `getReachAvoidSet()` respectively.
    - SReachTools/NAME: Succint description of the class/function
    - Horizontal separator
    - Detailed description
    - Usage example
    - If function, 
        - inputs 
            - If name-value pairs, list all possible names, values accepted, and
              optional/required
        - outputs
        - notes (Relevant papers)
    - If class, (make sure the components are linked, see
      [Matlab's help on help](https://www.mathworks.com/help/matlab/matlab_prog/create-help-for-classes.html))
        - properties 
            - List all except hidden and private
        - methods
            - List all except hidden, private, and overloaded MATLAB internals
            - List constructor
        - notes (Relevant papers)
    - Horizontal separator
    - Licensing information
    - 2 empty lines
- Subjective guidelines:
    - Try to align `=` when it improves readability.
- Follow goodpractices for commit message
    - Function/Class name changes must be specified in the commit message using
      a `->`.

### Glossary of acronyms used

To balance our pursuit of verbosity vs compactness, we have used throughout our
code base a list of acronyms. We list them below with their expansions:

- concat: concatenated 
- H: H\_matrix (see @LtiSystem/getConcatMats for notation)
- G: G\_matrix (see @LtiSystem/getConcatMats for notation)
- U: input space
- mat: matrix/matrices
- x0: initial state 
- cov: sigma
- mc: Monte-Carlo
- prob: probability
- thresh: threshold
- tol: tolerance
- n\_ : number of
- eff: effective
- aug: augmented
- dim: dimension
- underapprox: underapproximation/underapproximative
- Avoid the use of
    - Propositions like on, of, for
    - Terms like based, problem, object

## How to contribute?

- Make sure you have a GitHub account.
- Submit an issue, assuming one does not already exist.
    - Clearly describe the issue including steps to reproduce when it is a bug.
    - Make sure you fill in the earliest version that you know has the issue.
- Fork the repository on GitHub.

### Making changes

- Create a topic branch from where you want to base your work.
    - This is usually the master branch.
    - Only target release branches if you are certain your fix must be on that
      branch.
    - To quickly create a topic branch based on master, run 
      `git checkout -b fix/master/my_contribution master`. Please avoid working
      directly on the master branch.
- Follow the coding guidelines given above
- Make commits of logical units.
- Check for unnecessary whitespace with git diff --check before committing.
- Make sure your commit messages are in the proper format. If the commit
  addresses an issue filed, start the first line of the commit with the issue
  number in parentheses.
     ```
     (SR-1234) Make the example in CONTRIBUTING imperative and concrete
 
     Without this patch applied the example commit message in the CONTRIBUTING
     document is not a concrete example. This is a problem because the
     contributor is left to imagine what the commit message should look like
     based on a description rather than an example. This patch fixes the
     problem by making the example concrete and imperative.
 
     The first line is a real-life imperative statement with a ticket number
     from our issue tracker. The body describes the behavior without the patch,
     why this is a problem, and how the patch fixes the problem when applied.
     ``` 
- Make sure you have added the necessary tests for your changes.
- Run _all_ the tests to assure nothing else was accidentally broken.
- Use [erlang's otp commit
  guidelines](https://github.com/erlang/otp/wiki/writing-good-commit-messages)
  as a 
  reference to writing good commit messages.

### Submitting changes and revert policy

For now, send an email to [aby.vinod@gmail.com](aby.vinod@gmail.com) to coordinate this.

# Acknowledgements

This CONTRIBUTING.md file was derived from [Puppet's CONTRIBUTING.md](
ttps://github.com/puppetlabs/puppet/blob/master/CONTRIBUTING.md)
