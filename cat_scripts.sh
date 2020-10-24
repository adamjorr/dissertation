# USAGE: cat_scripts [-c //] [-s .cc] dir
# Concatenate scripts in the input directory.
# Each will have the comment characters appended to the top.

function usage(){
    echo "$0 [-c // -s .cc] dir"
}

CCHAR="//"
SUFFIX=".cc"

set -e
[ $# -eq 0 ] && usage && exit 1

while getopts "hc:s:" arg; do
    case $arg in
        c)
            CCHAR=${OPTARG}
            ;;
        s)
            SUFFIX=${OPTARG}
            ;;
        h | *)
            usage
            exit 0
    esac
done
shift $((OPTIND-1))

find $1 -name '*'${SUFFIX} -print0 | \
tac -s $'\0' | \
xargs -0 -n 1 -I{} echo cat "<(echo -e \\\n${CCHAR} FILE:{} \\\n\\\n)" {} | bash
