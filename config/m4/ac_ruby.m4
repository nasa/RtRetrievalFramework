# SYNOPSIS
#
#   AC_RUBY([version])
#
# DESCRIPTION
#
# if we find ruby on the system.

AC_DEFUN([AC_RUBY], 
[
AC_ARG_VAR([RUBY], [Override the Ruby executable to use,
           the default is just 'ruby'])
AC_ARG_VAR([RUBY_VERSION],[The installed Ruby
 	    version to use, for example '18'. This string
            will be appended to the Ruby interpreter
	    canonical name.])


if test "$RUBY" == ""; then
  RUBY=ruby
fi
AC_PATH_PROG([RUBY],[$RUBY[$RUBY_VERSION]],:,[$prefix/bin${PATH_SEPARATOR}$PATH])
if test "$ac_ruby_path" != ""; then
  RUBY=$ac_ruby_path
fi
])
