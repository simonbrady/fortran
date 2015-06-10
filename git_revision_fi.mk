# Tell make how to build git_revision.fi if it doesn't exist
git_revision.fi: git_revision
	
# Prevent unnecessary recompilation by only modifying git_revision.fi when the revision changes 
git_revision:
	git rev-parse HEAD | head -c 7 | xargs printf "character(len=*), parameter :: git_revision = '%s" > $@.tmp
ifeq ($(DEBUG), 1)
	echo -n "-debug" >> $@.tmp
endif
	echo "'" >> $@.tmp
	if cmp -s $@.tmp $@.fi; then \
		rm $@.tmp; \
	else \
		mv $@.tmp $@.fi; \
	fi

.PHONY: git_revision
