# This looks for the Lua and Luabind libraries. If we find them, we set 
# the Makefile conditional HAVE_LUA. We as set LUA_CPPFLAGS and LUA_LDFLAGS
# 
# To allow users to build there own copy of Lua, we also define
# BUILD_LUA

AC_DEFUN([AC_LUA],
[
AC_REQUIRE([MP_WITH_CURSES])
have_lua="no"
build_lua="no"
AC_ARG_WITH([lua],
        AS_HELP_STRING([--with-lua@<:@=DIR@:>@], [Build with Lua. You should have installed the thirdparty Lua libraries using the Lua target.]),
        [
    if test "$withval" = "no"; then
        want_lua="no"
    elif test "$withval" = "yes"; then
        want_lua="yes"
        ac_lua_path=""
    elif test "$withval" = "build"; then
        want_lua="yes"
        build_lua="yes"
        if test "x$prefix" = xNONE; then
           ac_lua_path="$ac_default_prefix"
        else
           ac_lua_path="$prefix"
        fi
    else
        want_lua="yes"
        ac_lua_path="$withval"
    fi
    ],
    [
    want_lua="yes"
    if test "$THIRDPARTY" = "build"; then
        want_lua="yes"
        build_lua="yes"
        if test "x$prefix" = xNONE; then
           ac_lua_path="$ac_default_prefix"
        else
           ac_lua_path="$prefix"
        fi
    fi
    ])

if test "x$want_lua" = "xyes"; then
        AC_MSG_CHECKING([for Lua library])
        succeeded=no
        if test "$ac_lua_path" != ""; then
            LUA_LDFLAGS="-L$ac_lua_path/lib -lluabind -llua -ldl"
            LUA_CPPFLAGS="-I$ac_lua_path/include"
            succeeded=yes
        else
            for ac_lua_path_tmp in $prefix $THIRDPARTY /groups/algorithm/commonlib/$commonlib_version/$compiler_name /usr /usr/local /opt /opt/local /sw ; do
                  if test -e "$ac_lua_path_tmp/include/luabind/luabind.hpp" && test -r "$ac_lua_path_tmp/include/luabind/luabind.hpp"; then
                      LUA_LDFLAGS="-L$ac_lua_path_tmp/lib -lluabind -llua -ldl"
                      LUA_CPPFLAGS="-I$ac_lua_path_tmp/include"
                      succeeded=yes
                      break;
                  fi
            done
        fi

        if test "$succeeded" != "yes" ; then
                AC_MSG_RESULT([no])
        else
                AC_MSG_RESULT([yes])
                AC_SUBST(LUA_CPPFLAGS)
                AC_SUBST(LUA_LDFLAGS)
                AC_DEFINE(HAVE_LUA,,[Defined if we have LUA])
                have_lua="yes"
        fi
fi
AM_CONDITIONAL([BUILD_LUA], [test "$build_lua" = "yes"])
AM_CONDITIONAL([HAVE_LUA], [test "$have_lua" = "yes"])
])