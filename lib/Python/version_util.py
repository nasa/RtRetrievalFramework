from __future__ import print_function
try:
    from builtins import str
except ImportError:
    # We might not have future installed yet, but this means we are using python 2.7 and 
    # can ignore this
    pass
import re
import os
import sys
import subprocess

class VersionError(Exception):
    pass

def _which(program):
    import os
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def _get_svn_version(source_path):
    "Get version information for a Subversion checkout. Due to the slow speed of svnversion, this will not return any indication if the source has been modified."

    # Make sure we can query commandline tools
    if _which("svn") == None:
        raise VersionError("svn command not found")

    svn_path = None
    svn_revision = None
    try:
        info_process = subprocess.Popen(["svn", "info", source_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        info_process.wait()

        for ver_line in info_process.stdout.readlines():
            if ver_line.find(b'URL:') >= 0:
                svn_path = ver_line.replace(b'URL: ',b'').strip()
            if ver_line.find(b'Revision:') >= 0:
                svn_revision = ver_line.replace(b'Revision: ',b'').strip()
        info_process.stdout.close()
    except Exception as exc:
        raise VersionError("Error querying svn command: %s" % str(exc))

    if svn_path == None:
        return None

    tag_match = re.search(b'tags/([^/]+)(/[^/]*)*$', svn_path)
    if tag_match:
        tag_name = tag_match.group(1)
    else:
        tag_name = None

    branch_match = re.search(b'branches/([^/]+)(/[^/]*)*$', svn_path)
    if branch_match:
        branch_name = branch_match.group(1)
    else:
        branch_name = None

    tag_branch = ''
    if tag_name != None:
        tag_branch = tag_name
    elif branch_name != None:
        tag_branch = branch_name

    rev_str = None
    if len(tag_branch) > 0:
        rev_str = tag_branch + b'-' + svn_revision
    else:
        rev_str = svn_revision

    return b'SVN-' + rev_str

def _get_git_version(source_path):
    "Get version information for a git clone"

    # Make sure we can query commandline tools
    if _which("git") == None:
        raise VersionError("git command not found")  
    
    prev_dir = os.getcwd()
    # git rev-parse works from cwd
    os.chdir(source_path) 
    try:
        git_process = subprocess.Popen(["git", "rev-parse", "--verify", "HEAD", "--short"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        git_process.wait()

        git_ver = git_process.stdout.readline().strip()
    except OSError:
        git_ver = None

    # Try and get branch version
    git_branch = b'(no branch)'
    try:
        git_process = subprocess.Popen(["git", "symbolic-ref", "HEAD"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        git_process.wait()
        
        if git_process.returncode == 0:
            git_branch = git_process.stdout.readline().strip().replace(b"refs/heads/", b"")
    except OSError:
        git_branch = None

    # See if working directorty has been modified
    wdir_modified = b'?'
    try:
        git_process = subprocess.Popen(["git", "status", "--porcelain", "-uno"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        git_process.wait()

        status_ret = git_process.stdout.readline().strip()
        if len(status_ret) > 0:
            wdir_modified = b'M'
        else:
            wdir_modified = b''
    except OSError:
        wdir_modified = None

    os.chdir(prev_dir)
    return b"GIT-" + git_branch + b'-' + git_ver + wdir_modified

def source_version(source_path):
    "Retrieves version information for a GIT or Subversion directory. The returned string includes branch names as well as a revision number or hash. We first check if the path is a git clone, failing that we try subversion."

    if not os.path.exists(source_path):
        return None

    # Detect if we are using Git for referenced file
    prev_dir = os.getcwd()
    os.chdir(source_path)
    if _which("git") != None:
        # Send stdout to a PIPE so it is not seen
        ret_code = subprocess.call(["git", "rev-parse", "--git-dir"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    else:
        ret_code = None
        
    os.chdir(prev_dir)
    
    if ret_code != None and ret_code == 0:
        return _get_git_version(source_path)
    else:
        return _get_svn_version(source_path)

def binary_version(binary_filename):
    "Query a L2 binary to get the compiled in version information"

    process = subprocess.Popen([binary_filename, "-v"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    ret_code = process.wait()

    if ret_code == 0:
        # L2 spits out several lines of version information when called with -v
        # If the first line is not as expected then this is not a l2 binary
        first_line = process.stdout.readline()
        if not first_line.find(b"Major version:") == 0:
            return None

        major_version = first_line.split(b":")[1].strip()
        cm_version = process.stdout.readline().split(b":")[1].strip()
        lua_version = process.stdout.readline()
        if len(lua_version) > 0:
            lua_version = lua_version.split(b":")[1].strip()
        else:
            lua_version = None

        return major_version, cm_version, lua_version
    else:
        return None

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please supply a path to a Git or Subversion repository", file=sys.stderr)
        sys.exit(1)

    if sys.version_info > (3,):
        print(source_version(sys.argv[1]).decode('utf-8'))
    else:
        print(source_version(sys.argv[1]))
