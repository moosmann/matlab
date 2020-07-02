function com_id =  git_commit_id()
% Return the commit ID of the currently used git repository.

cur_dir = pwd;
cd( userpath );
[~, com_id] = unix( 'git rev-parse HEAD' );
com_id(end) = [];
cd( cur_dir );