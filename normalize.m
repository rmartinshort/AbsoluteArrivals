DO FILE WILD $1
    read $FILE
    lh depmin depmax
    setbb a (max &1,depmax)
    setbb b (min &1,depmin)
    div (max %a (abs %b))
    write over
ENDDO
