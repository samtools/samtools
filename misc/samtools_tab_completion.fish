################################################################################
#
# The MIT License
#
# Copyright 2025 LunarEclipse <luna@lunareclipse.zone>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the “Software”), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
################################################################################

# Collect subcommands and their descriptions into a variable for later use
# Exclude "help" because "samtools help help" won't work (and wouldn't parse correctly)
set -l __fish_samtools_subcommands (samtools 2>&1 \
    | string match --entire --regex '^   ' \
    | string trim \
    | string replace --regex '\s{2,}' \t \
    | string match --invert 'help*' \
    | string escape \
)
# Add "help" back to an extra helper variable
set -l __fish_samtools_subcommands_with_help $__fish_samtools_subcommands
set -la __fish_samtools_subcommands_with_help "help\\t'print help for a subcommand'"

# Subcommands
complete -c samtools --no-files --condition __fish_use_subcommand -a "$__fish_samtools_subcommands_with_help"
complete -c samtools --no-files --condition '__fish_seen_subcommand_from help' --require-parameter -a "$__fish_samtools_subcommands"

# Generation from manpages doesn't work well for root options
complete -c samtools --condition __fish_use_subcommand -l help -d 'Print help for a subcommand.' --no-files --require-parameter -a "$__fish_samtools_subcommands"
complete -c samtools --condition __fish_use_subcommand -l version -d 'Display the version numbers and copyright information for samtools and the libraries it uses.'
complete -c samtools --condition __fish_use_subcommand -l version-only -d 'Display the full version number in a machine-readable format.'

# Determine if and where the samtools manpages are installed
which python3 &>/dev/null
and python3 -c "import sys, os; sys.path.append('$__fish_data_dir/tools'); from create_manpage_completions import get_paths_from_man_locations; print(os.path.dirname([x for x in get_paths_from_man_locations() if 'samtools.1' in x][0]))" 2>/dev/null | read -l __fish_samtools_man_path
and set -l __fish_samtools_manpages "$__fish_samtools_man_path"/samtools-*

# Generate completions for subcommands from manpages if they're available
# Otherwise fall back to the approach bash completions take
if count $__fish_samtools_manpages &>/dev/null
    "$__fish_data_dir"/tools/create_manpage_completions.py --stdout $__fish_samtools_manpages \
        | string replace --regex -- '-c samtools-(\w+)' '-c samtools -n \'__fish_seen_subcommand_from $1\'' \
        | read -lz __fish_samtools_subcommand_completions
    eval "$__fish_samtools_subcommand_completions"
else
    for subcommand in (string unescape $__fish_samtools_subcommands | string match --regex '.*\t' | string trim);
        for line in $(samtools help "$subcommand" 2>&1);
            # Note: this will only use the last version of a long argument
            # Note: this will ignore long commands containing braces or other special characters
            # https://regex101.com/r/4IYuNh/5
            string match --quiet --regex '(?J)(?:(?:[ ,]-(?<short>[\w@?])[ ,])*(?=.*?)(?:[ ,]--(?<long>[-\w_]+)(?:[ ,]|$))|(?:^|Options?:?) *(?:[ ,]?(?:\(or )?-(?<short>[\w@?]{1,2})[ ,)])+)' "$line"; or continue

            if count $short &>/dev/null
                set arg_short --short-option "$short"
                test (string length $short) -gt 1; and set arg_short --old-option "$short"
            end
            count $long &>/dev/null; and set arg_long --long-option "$long"

            complete -c samtools --condition "__fish_seen_subcommand_from $subcommand" $arg_short $arg_long
        end
    end
end

