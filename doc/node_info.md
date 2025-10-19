[glucarel@login07 HPC-leonardo]$ srun --partition=dcgp_usr_prod \
>      --account=uTS25_Tornator_0 \
>      --nodes=1 \
>      --ntasks=1 \
>      --cpus-per-task=112 \
>      --time=00:15:00 \
>      --exclusive \
>      --pty bash
srun: no gres/tmpfs specified, using default: gres/tmpfs:10g
srun: job 21697375 queued and waiting for resources
srun: job 21697375 has been allocated resources
[glucarel@lrdn4293 HPC-leonardo]$ numactl --hardware
available: 8 nodes (0-7)
node 0 cpus: 0 1 2 3 4 5 6 7 8 9 10 11 12 13
node 0 size: 63565 MB
node 0 free: 61356 MB
node 1 cpus: 14 15 16 17 18 19 20 21 22 23 24 25 26 27
node 1 size: 64508 MB
node 1 free: 64112 MB
node 2 cpus: 28 29 30 31 32 33 34 35 36 37 38 39 40 41
node 2 size: 64508 MB
node 2 free: 62688 MB
node 3 cpus: 42 43 44 45 46 47 48 49 50 51 52 53 54 55
node 3 size: 64508 MB
node 3 free: 64015 MB
node 4 cpus: 56 57 58 59 60 61 62 63 64 65 66 67 68 69
node 4 size: 64464 MB
node 4 free: 62968 MB
node 5 cpus: 70 71 72 73 74 75 76 77 78 79 80 81 82 83
node 5 size: 64508 MB
node 5 free: 64133 MB
node 6 cpus: 84 85 86 87 88 89 90 91 92 93 94 95 96 97
node 6 size: 64508 MB
node 6 free: 64076 MB
node 7 cpus: 98 99 100 101 102 103 104 105 106 107 108 109 110 111
node 7 size: 64505 MB
node 7 free: 63087 MB
node distances:
node   0   1   2   3   4   5   6   7 
  0:  10  12  12  12  21  21  21  21 
  1:  12  10  12  12  21  21  21  21 
  2:  12  12  10  12  21  21  21  21 
  3:  12  12  12  10  21  21  21  21 
  4:  21  21  21  21  10  12  12  12 
  5:  21  21  21  21  12  10  12  12 
  6:  21  21  21  21  12  12  10  12 
  7:  21  21  21  21  12  12  12  10 
[glucarel@lrdn4293 HPC-leonardo]$ lscpu -e=CPU,CORE,SOCKET,NODE,ONLINE
CPU CORE SOCKET NODE ONLINE
0   0    0      0    yes
1   1    0      0    yes
2   2    0      0    yes
3   3    0      0    yes
4   4    0      0    yes
5   5    0      0    yes
6   6    0      0    yes
7   7    0      0    yes
8   8    0      0    yes
9   9    0      0    yes
10  10   0      0    yes
11  11   0      0    yes
12  12   0      0    yes
13  13   0      0    yes
14  14   0      1    yes
15  15   0      1    yes
16  16   0      1    yes
17  17   0      1    yes
18  18   0      1    yes
19  19   0      1    yes
20  20   0      1    yes
21  21   0      1    yes
22  22   0      1    yes
23  23   0      1    yes
24  24   0      1    yes
25  25   0      1    yes
26  26   0      1    yes
27  27   0      1    yes
28  28   0      2    yes
29  29   0      2    yes
30  30   0      2    yes
31  31   0      2    yes
32  32   0      2    yes
33  33   0      2    yes
34  34   0      2    yes
35  35   0      2    yes
36  36   0      2    yes
37  37   0      2    yes
38  38   0      2    yes
39  39   0      2    yes
40  40   0      2    yes
41  41   0      2    yes
42  42   0      3    yes
43  43   0      3    yes
44  44   0      3    yes
45  45   0      3    yes
46  46   0      3    yes
47  47   0      3    yes
48  48   0      3    yes
49  49   0      3    yes
50  50   0      3    yes
51  51   0      3    yes
52  52   0      3    yes
53  53   0      3    yes
54  54   0      3    yes
55  55   0      3    yes
56  56   1      4    yes
57  57   1      4    yes
58  58   1      4    yes
59  59   1      4    yes
60  60   1      4    yes
61  61   1      4    yes
62  62   1      4    yes
63  63   1      4    yes
64  64   1      4    yes
65  65   1      4    yes
66  66   1      4    yes
67  67   1      4    yes
68  68   1      4    yes
69  69   1      4    yes
70  70   1      5    yes
71  71   1      5    yes
72  72   1      5    yes
73  73   1      5    yes
74  74   1      5    yes
75  75   1      5    yes
76  76   1      5    yes
77  77   1      5    yes
78  78   1      5    yes
79  79   1      5    yes
80  80   1      5    yes
81  81   1      5    yes
82  82   1      5    yes
83  83   1      5    yes
84  84   1      6    yes
85  85   1      6    yes
86  86   1      6    yes
87  87   1      6    yes
88  88   1      6    yes
89  89   1      6    yes
90  90   1      6    yes
91  91   1      6    yes
92  92   1      6    yes
93  93   1      6    yes
94  94   1      6    yes
95  95   1      6    yes
96  96   1      6    yes
97  97   1      6    yes
98  98   1      7    yes
99  99   1      7    yes
100 100  1      7    yes
101 101  1      7    yes
102 102  1      7    yes
103 103  1      7    yes
104 104  1      7    yes
105 105  1      7    yes
106 106  1      7    yes
107 107  1      7    yes
108 108  1      7    yes
109 109  1      7    yes
110 110  1      7    yes
111 111  1      7    yes
[glucarel@lrdn4293 HPC-leonardo]$ lscpu | egrep 'L1d|L1i|L2|L3'
L1d cache:           48K
L1i cache:           32K
L2 cache:            2048K
L3 cache:            107520K
[glucarel@lrdn4293 HPC-leonardo]$ for i in /sys/devices/system/cpu/cpu0/cache/index*; do
>   echo -n "$i: "
>   cat $i/level $i/type $i/size
> done
/sys/devices/system/cpu/cpu0/cache/index0: 1
Data
48K
/sys/devices/system/cpu/cpu0/cache/index1: 1
Instruction
32K
/sys/devices/system/cpu/cpu0/cache/index2: 2
Unified
2048K
/sys/devices/system/cpu/cpu0/cache/index3: 3
Unified
107520K