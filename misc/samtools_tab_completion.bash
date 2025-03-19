################################################################################
# Copyright (c) 2016 Genome Research Ltd. 
# 
# Author: George Hall <gh10@sanger.ac.uk> 
# 
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
################################################################################

# This script determines whether the user is trying to complete a subcommand or 
# a long-option, then determines and displays possible completions. For commands
# with long-options, the initial '-' must be first be typed in order to trigger 
# the completion mechanism. 

# This file must be sourced for the tab completions to work - it is advisable to 
# source it from your .bashrc. By default, it will only perform tab completion for
# a binary called 'samtools'. To enable tab completion for other versions of
# Samtools, or for Samtools binaries with different names, simply append 
# 'complete -F _samtools_options <name of binary>' to the end of this file. This 
# file will then need to be sourced again. 

_samtools_options() 
{

	SAMTOOLS_BIN="${COMP_WORDS[0]}"
	local CUR CURRENT_SUB OPTS
	COMPREPLY=()
	CUR="${COMP_WORDS[COMP_CWORD]}"

	if [[ $COMP_CWORD -eq 1 ]] ; then

		# If on the first word, generate possible subcommands, and tab complete those

		OPTS=$($SAMTOOLS_BIN 2>&1 | grep '^   ' | awk '{print $1}')

		COMPREPLY=($(compgen -W "${OPTS}" -- ${CUR}))
		return 0

	else 

		# Complete long-options (if available) for current subcommand

		CURRENT_SUB="${COMP_WORDS[1]}"

		OPTS=$($SAMTOOLS_BIN $CURRENT_SUB 2>&1 | grep -oh "\(\-\-\)\([[:alnum:]]\|\-\)* " | \
			xargs -I @ printf -- @)
		
		if [[ ${CUR} == -* ]] ; then
			COMPREPLY=($(compgen -W "${OPTS}" -- ${CUR}))
			return 0
		else
			# Assume the user wants normal file name completion
			COMPREPLY=($(compgen -o default))
		fi

	fi

}

complete -F _samtools_options samtools

