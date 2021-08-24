#
# Python module for TreeFix-TP
#

#=============================================================================
# constants

PROGRAM_NAME = "TreeFix-TP"
PROGRAM_VERSION_MAJOR = 1
PROGRAM_VERSION_MINOR = 2
PROGRAM_VERSION_RELEASE = 2
PROGRAM_VERSION = (PROGRAM_VERSION_MAJOR,
                   PROGRAM_VERSION_MINOR,
                   PROGRAM_VERSION_RELEASE)

if PROGRAM_VERSION_RELEASE != 0:
    PROGRAM_VERSION_TEXT = "%d.%d.%d" % (PROGRAM_VERSION_MAJOR,
                                         PROGRAM_VERSION_MINOR,
                                         PROGRAM_VERSION_RELEASE)
else:
    PROGRAM_VERSION_TEXT = "%d.%d" % (PROGRAM_VERSION_MAJOR,
                                      PROGRAM_VERSION_MINOR)
