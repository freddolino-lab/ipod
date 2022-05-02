BEGIN {
    OFS=""
    FS="\t"
}
{
    p2 = NR+2
    # if line number plus 2 is divisible by 4, it's the seq
    if (p2 % 4 == 0) {
        print substr($1,9), $2
    # if line number is divisible by 4, it's the qual
    } else if (NR % 4 == 0) {
        print substr($1,9), $2
    # any other line is either name or "+"
    } else {
        print $2
    }
}

