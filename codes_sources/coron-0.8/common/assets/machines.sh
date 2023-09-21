# This file is sourced into every launch script.
# Put your machine info in the list in alphabetical order.

if [ $AUTO_MEM = "1" ]; then
   case "$HOSTNAME" in
      "BARRACUDA"            )   JAVA_MAX_HEAP_SIZE=1150m;;
      "coutures"             )   JAVA_MAX_HEAP_SIZE=512m;;
      "debrecen"             )   JAVA_MAX_HEAP_SIZE=1024m;;
      "hagrid.loria.fr"      )   JAVA_MAX_HEAP_SIZE=1024m;;
      "hoth"                 )   JAVA_MAX_HEAP_SIZE=2600m;;
      "kraken"	    	        )   JAVA_MAX_HEAP_SIZE=2660m;;
      "servbioinfo.loria.fr" )   JAVA_MAX_HEAP_SIZE=2630m;;
      *                      )   JAVA_MAX_HEAP_SIZE=DEFAULT;;
   esac
else
   JAVA_MAX_HEAP_SIZE=DEFAULT
fi
