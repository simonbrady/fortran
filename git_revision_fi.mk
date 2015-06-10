# Tell make how to build git_revision.fi if it doesn't exist
git_revision.fi: git_revision
	
# Prevent unnecessary recompilation by only modifying git_revision.fi when the revision changes 
git_revision:
	rev=$$(git rev-parse --short HEAD 2>/dev/null || echo "unknown"); \
	printf "character(len=*), parameter :: git_revision = '%s'\n" $${rev} > $@.tmp; \
	if cmp -s $@.tmp $@.fi; then \
		rm $@.tmp; \
	else \
		mv $@.tmp $@.fi; \
	fi

.PHONY: git_revision
