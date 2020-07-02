function d = git_commit_date( commit_ID )
% Return the date of a commit. Default: use last commit.


%% NOT WORKING
if nargin < 1
    commit_ID = git_commit_id();
end

s = sprintf( 'cd %s; git show -s --format=%%ci %s', userpath, commit_ID );
disp(s)
%s = sprintf( 'cd %s; pwd', userpath);
unix( s );
