#!/bin/sh
set -e
java -Xms48G -Xmx48G -jar `dirname $0`/juicebox_tools.jar $*
