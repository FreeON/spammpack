#!/bin/bash
#
# This script updates the file git_commit_tag and returns 1 if it updated the
# file using git and 0 otherwise.

GIT_COMMIT_FILE=git_commit_tag

return_value=0

git --version > /dev/null 2>&1

if test $? -eq 0; then
	BRANCH=`git branch 2>&1 | egrep '^\*' | sed -e 's/^\* //' | sed -e 's/[()]//g' | sed -e 's/ /_/g'`
	COMMIT=`git show --pretty=oneline HEAD | head -n1 | awk '{print $1}'`
	echo "branch=$BRANCH" > ${GIT_COMMIT_FILE}.temp
	echo "git_commit_tag=$COMMIT" >> ${GIT_COMMIT_FILE}.temp
	diff --brief ${GIT_COMMIT_FILE} ${GIT_COMMIT_FILE}.temp > /dev/null 2>&1
	if test $? -eq 0; then
		echo "not updating ${GIT_COMMIT_FILE}, no change"
	else
		echo "updating ${GIT_COMMIT_FILE}"
		mv -f ${GIT_COMMIT_FILE}.temp ${GIT_COMMIT_FILE}
		return_value=1
	fi
else
  if test -f ${GIT_COMMIT_FILE}; then
    echo "no git, hence ${GIT_COMMIT_FILE} is assumed to be correct"
  else
    echo "no ${GIT_COMMIT_FILE}, creating one"
    echo "branch=unknown" > ${GIT_COMMIT_FILE}
    echo "commit=unknown" >> ${GIT_COMMIT_FILE}
  fi
fi

rm -f ${GIT_COMMIT_FILE}.temp

exit $return_value
