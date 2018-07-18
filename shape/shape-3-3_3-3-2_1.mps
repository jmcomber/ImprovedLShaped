* Problem:
* Class:      MIP
* Rows:       6
* Columns:    9 (9 integer, 9 binary)
* Non-zeros:  36
* Format:     Fixed MPS
*
NAME
ROWS
 N  R0000000
 G  c1
 G  c2
 G  c3
 G  c4
 G  c5
 G  c6
COLUMNS
    M0000001  'MARKER'                 'INTORG'
    p_00      R0000000             5   c1                  22
    p_00      c2                  60   c3                   2
    p_01      R0000000            46   c1                   1
    p_01      c2                  40   c3                  12
    p_02      R0000000            84   c1                  58
    p_02      c2                  90   c3                  89
    x_00      R0000000            34   c1                  94
    x_00      c2                  70   c3                  12
    x_00      c4                  32   c5                  32
    x_00      c6                  38
    x_01      R0000000            57   c1                  35
    x_01      c2                  23   c3                  76
    x_01      c4                  44   c5                  99
    x_01      c6                  57
    x_02      R0000000             1   c1                  90
    x_02      c2                  97   c3                  30
    x_02      c4                  67   c5                  68
    x_02      c6                  48
    q_00      R0000000             1   c4                  65
    q_00      c5                  71   c6                  16
    q_01      R0000000             1   c4                  43
    q_01      c5                  78   c6                  88
    q_02      R0000000             1   c4                  51
    q_02      c5                  82   c6                  70
    M0000002  'MARKER'                 'INTEND'
RHS
    RHS1      c1                 225   c2                 285
    RHS1      c3                 165   c4                 226
    RHS1      c5                 322   c6                 237
BOUNDS
 UP BND1      p_00                 1
 UP BND1      p_01                 1
 UP BND1      p_02                 1
 UP BND1      x_00                 1
 UP BND1      x_01                 1
 UP BND1      x_02                 1
 UP BND1      q_00                 1
 UP BND1      q_01                 1
 UP BND1      q_02                 1
ENDATA
