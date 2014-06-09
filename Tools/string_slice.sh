function string_slice {
    STRING="$1"
    declare -i LENGTH="${#STRING}"
    declare -i START="$2"
    declare -i END="$3"
    if [ $START -lt 0 ]; then
        START=$[ $LENGTH + $START ]
    fi
    if [ $END -le 0 ]; then
        END=$[ $LENGTH + $END ]
    fi
    START=$[ $START + 1 ]
    (echo "$STRING" | cut -c $START-$END) 2> /dev/null
}
