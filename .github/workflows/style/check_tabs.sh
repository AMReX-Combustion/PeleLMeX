#!/usr/bin/env bash

set -eu -o pipefail

find . -type d \( -name .git \
                  -o -name Submodules \
               \) -prune -o \
       -type f \( \( -name "*.H" -o -name "*.h" -o -name "*.hh" -o -name "*.hpp" \
                  -o -name "*.c" -o -name "*.cc" -o -name "*.cpp" -o -name "*.cxx" \
                  -o -name "*.py" \
                  -o -name "*.md" -o -name "*.rst" \
                  -o -name "*.sh" \
                  -o -name "*.tex" \
                  -o -name "*.txt" \
                  -o -name "*.yml" \) \
                 -a \( ! -name "*.tab.h" -a ! -name "*.tab.cpp" \
                    -a ! -name "*.lex.h" -a ! -name "*.lex.cpp" \) \
               \) \
    -exec grep -Iq . {} \; \
    -exec sed -i 's/\t/\ \ \ \ /g' {} +

gitdiff=`git diff`

if [ -z "$gitdiff" ]
then
    exit 0
else
    echo -e "\nTabs are not allowed. Changes suggested by"
    echo -e "  ${0}\n"
    git --no-pager diff
    echo ""
    exit 1
fi
