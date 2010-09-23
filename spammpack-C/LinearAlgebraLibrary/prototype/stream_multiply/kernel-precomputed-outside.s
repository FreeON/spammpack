GNU gdb (Gentoo 7.1 p1) 7.1
Copyright (C) 2010 Free Software Foundation, Inc.
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.  Type "show copying"
and "show warranty" for details.
This GDB was configured as "x86_64-pc-linux-gnu".
For bug reporting instructions, please see:
<http://bugs.gentoo.org/>.
Dump of assembler code for function stream_kernel:
74	{
   0x0000000000000000 <+0>:	sub    $0x1e8,%rsp

75	  short int i;
76	  unsigned int stream_index;
77	  unsigned int max_stream_index;
78
79	  float *restrict A;
80	  float *restrict B;
81	  float *restrict C;
82
83	  float *restrict norm;
84
85	  __m128 alpha_row;
86
87	  __m128 A_element;
88	  __m128 B_row;
89	  __m128 C_row[4];
90
91	  char norm_product[16][16];
92
93	  /* Divide number of stream elements by 64 to simulate stride of 64. */
94	  max_stream_index = number_stream_elements/64;
   0x000000000000000d <+13>:	shr    $0x6,%edi

95
96	  alpha_row = _mm_set1_ps(alpha);
   0x0000000000000007 <+7>:	shufps $0x0,%xmm0,%xmm0

97
98	  for (stream_index = 0; stream_index < max_stream_index; stream_index++)
   0x000000000000000b <+11>:	xor    %edx,%edx
   0x0000000000000010 <+16>:	test   %edi,%edi
   0x0000000000000012 <+18>:	jbe    0x5bda <stream_kernel+23514>
   0x0000000000000018 <+24>:	xor    %eax,%eax
   0x000000000000001a <+26>:	movss  %xmm1,0x1c0(%rsp)
   0x0000000000005b89 <+23433>:	inc    %eax
   0x0000000000005b8b <+23435>:	mov    %eax,%edx
   0x0000000000005b95 <+23445>:	mov    %edx,%eax
   0x0000000000005b97 <+23447>:	cmp    %edi,%eax
   0x0000000000005bd4 <+23508>:	jb     0x1a <stream_kernel+26>

99	  {
100	    /* Load pointers to matrix data blocks. */
101	    A = multiply_stream[stream_index].A_block;
   0x0000000000000023 <+35>:	imul   $0x98,%rdx,%rdx
   0x000000000000002a <+42>:	mov    (%rsi,%rdx,1),%r8

102	    B = multiply_stream[stream_index].B_block;
   0x000000000000002e <+46>:	mov    0x8(%rsi,%rdx,1),%rcx

103	    C = multiply_stream[stream_index].C_block;
   0x0000000000000033 <+51>:	mov    0x10(%rsi,%rdx,1),%r9

104	    norm = multiply_stream[stream_index].norm;
105
106	    /* Calculate norms. */
107	    norm_product[0][0] = (norm[0]*norm[16] >= tolerance);
   0x0000000000000038 <+56>:	movss  0x18(%rdx,%rsi,1),%xmm4
   0x000000000000003e <+62>:	movss  0x58(%rdx,%rsi,1),%xmm1
   0x000000000000007a <+122>:	movaps %xmm4,%xmm7
   0x000000000000007d <+125>:	mulss  %xmm1,%xmm7
   0x000000000000008c <+140>:	movss  %xmm7,0x1c8(%rsp)

108	    norm_product[1][4] = (norm[1]*norm[20] >= tolerance);
   0x0000000000000044 <+68>:	movss  0x1c(%rdx,%rsi,1),%xmm3
   0x0000000000000057 <+87>:	movss  0x68(%rdx,%rsi,1),%xmm13
   0x0000000000000081 <+129>:	movaps %xmm3,%xmm6
   0x0000000000000087 <+135>:	mulss  %xmm13,%xmm6
   0x00000000000000a6 <+166>:	movss  %xmm6,0x1b0(%rsp)

109	    norm_product[2][8] = (norm[2]*norm[24] >= tolerance);
   0x000000000000004a <+74>:	movss  0x20(%rdx,%rsi,1),%xmm2
   0x000000000000005e <+94>:	movss  0x78(%rdx,%rsi,1),%xmm11
   0x0000000000000084 <+132>:	movaps %xmm2,%xmm5
   0x00000000000000a1 <+161>:	mulss  %xmm11,%xmm5
   0x00000000000000be <+190>:	movss  %xmm5,0x1a8(%rsp)

110	    norm_product[3][12] = (norm[3]*norm[28] >= tolerance);
   0x0000000000000050 <+80>:	movss  0x24(%rdx,%rsi,1),%xmm10
   0x0000000000000065 <+101>:	movss  0x88(%rdx,%rsi,1),%xmm12
   0x0000000000000095 <+149>:	movaps %xmm10,%xmm7
   0x00000000000000b9 <+185>:	mulss  %xmm12,%xmm7
   0x00000000000000cf <+207>:	movss  %xmm7,0x1a0(%rsp)

111	    norm_product[0][1] = (norm[0]*norm[17] >= tolerance);
   0x00000000000000af <+175>:	movaps %xmm4,%xmm6
   0x00000000000000d8 <+216>:	movss  0x5c(%rdx,%rsi,1),%xmm7
   0x00000000000000de <+222>:	mulss  %xmm7,%xmm6
   0x00000000000000f8 <+248>:	movss  %xmm6,0x1b8(%rsp)

112	    norm_product[1][5] = (norm[1]*norm[21] >= tolerance);
   0x000000000000006f <+111>:	movss  0x6c(%rdx,%rsi,1),%xmm15
   0x0000000000000076 <+118>:	movaps %xmm0,(%rsp)
   0x0000000000000099 <+153>:	movaps %xmm3,%xmm9
   0x00000000000000ca <+202>:	mulss  %xmm15,%xmm9
   0x00000000000000e2 <+226>:	movss  %xmm9,0x190(%rsp)

113	    norm_product[2][9] = (norm[2]*norm[25] >= tolerance);
   0x000000000000009d <+157>:	movaps %xmm2,%xmm14
   0x00000000000000ec <+236>:	movss  0x7c(%rdx,%rsi,1),%xmm9
   0x00000000000000f3 <+243>:	mulss  %xmm9,%xmm14
   0x000000000000010b <+267>:	movss  %xmm14,0x180(%rsp)

114	    norm_product[3][13] = (norm[3]*norm[29] >= tolerance);
   0x00000000000000b2 <+178>:	movaps %xmm10,%xmm8
   0x0000000000000115 <+277>:	movss  0x8c(%rdx,%rsi,1),%xmm14
   0x000000000000011f <+287>:	mulss  %xmm14,%xmm8
   0x0000000000000136 <+310>:	movss  %xmm8,0x178(%rsp)

115	    norm_product[0][2] = (norm[0]*norm[18] >= tolerance);
   0x00000000000000b6 <+182>:	movaps %xmm4,%xmm0
   0x0000000000000101 <+257>:	movss  0x60(%rdx,%rsi,1),%xmm6
   0x0000000000000107 <+263>:	mulss  %xmm6,%xmm0
   0x0000000000000124 <+292>:	movss  %xmm0,0x198(%rsp)

116	    norm_product[1][6] = (norm[1]*norm[22] >= tolerance);
   0x00000000000000c7 <+199>:	movaps %xmm3,%xmm5
   0x0000000000000140 <+320>:	movss  0x70(%rdx,%rsi,1),%xmm8
   0x0000000000000147 <+327>:	mulss  %xmm8,%xmm5
   0x000000000000014c <+332>:	movss  %xmm5,0x168(%rsp)

117	    norm_product[2][10] = (norm[2]*norm[26] >= tolerance);
   0x000000000000012d <+301>:	movss  0x80(%rdx,%rsi,1),%xmm0
   0x0000000000000155 <+341>:	movaps %xmm2,%xmm5
   0x0000000000000158 <+344>:	mulss  %xmm0,%xmm5
   0x000000000000015c <+348>:	movss  %xmm5,0x158(%rsp)

118	    norm_product[3][14] = (norm[3]*norm[30] >= tolerance);
   0x0000000000000165 <+357>:	movaps %xmm10,%xmm5
   0x0000000000000169 <+361>:	mulss  0x90(%rdx,%rsi,1),%xmm5
   0x0000000000000172 <+370>:	movss  %xmm5,0x150(%rsp)

119	    norm_product[0][3] = (norm[0]*norm[19] >= tolerance);
   0x000000000000017b <+379>:	movss  0x64(%rdx,%rsi,1),%xmm5
   0x0000000000000181 <+385>:	mulss  %xmm5,%xmm4
   0x0000000000000185 <+389>:	movss  %xmm4,0x1d0(%rsp)

120	    norm_product[1][7] = (norm[1]*norm[23] >= tolerance);
   0x000000000000018e <+398>:	movss  0x74(%rdx,%rsi,1),%xmm4
   0x0000000000000194 <+404>:	mulss  %xmm4,%xmm3
   0x0000000000000198 <+408>:	movss  %xmm3,0x140(%rsp)

121	    norm_product[2][11] = (norm[2]*norm[27] >= tolerance);
   0x00000000000001a1 <+417>:	movss  0x84(%rdx,%rsi,1),%xmm3
   0x00000000000001aa <+426>:	mulss  %xmm3,%xmm2
   0x00000000000001ae <+430>:	movss  %xmm2,0x138(%rsp)

122	    norm_product[3][15] = (norm[3]*norm[31] >= tolerance);
   0x00000000000001b7 <+439>:	movss  0x94(%rdx,%rsi,1),%xmm2
   0x00000000000001c0 <+448>:	mulss  %xmm2,%xmm10
   0x00000000000001c5 <+453>:	movss  %xmm10,0x130(%rsp)

123	    norm_product[4][0] = (norm[4]*norm[16] >= tolerance);
   0x00000000000001cf <+463>:	movss  0x28(%rdx,%rsi,1),%xmm10
   0x00000000000001d6 <+470>:	mulss  %xmm10,%xmm1
   0x00000000000001db <+475>:	movss  %xmm1,0x170(%rsp)

124	    norm_product[5][4] = (norm[5]*norm[20] >= tolerance);
   0x00000000000001e4 <+484>:	movss  0x2c(%rdx,%rsi,1),%xmm1
   0x00000000000001ea <+490>:	mulss  %xmm1,%xmm13
   0x00000000000001ef <+495>:	movss  %xmm13,0x120(%rsp)

125	    norm_product[6][8] = (norm[6]*norm[24] >= tolerance);
   0x00000000000001f9 <+505>:	movaps %xmm11,%xmm13
   0x00000000000001fd <+509>:	movss  0x30(%rdx,%rsi,1),%xmm11
   0x0000000000000204 <+516>:	mulss  %xmm11,%xmm13
   0x0000000000000209 <+521>:	movss  %xmm13,0x118(%rsp)

126	    norm_product[7][12] = (norm[7]*norm[28] >= tolerance);
   0x0000000000000213 <+531>:	movss  0x34(%rdx,%rsi,1),%xmm13
   0x000000000000021a <+538>:	mulss  %xmm13,%xmm12
   0x000000000000021f <+543>:	movss  %xmm12,0x110(%rsp)

127	    norm_product[4][1] = (norm[4]*norm[17] >= tolerance);
   0x0000000000000229 <+553>:	movaps %xmm7,%xmm12
   0x000000000000022d <+557>:	mulss  %xmm10,%xmm12
   0x0000000000000232 <+562>:	movss  %xmm12,0x160(%rsp)

128	    norm_product[5][5] = (norm[5]*norm[21] >= tolerance);
   0x000000000000023c <+572>:	movaps %xmm15,%xmm12
   0x0000000000000240 <+576>:	mulss  %xmm1,%xmm12
   0x0000000000000245 <+581>:	movss  %xmm12,0x100(%rsp)

129	    norm_product[6][9] = (norm[6]*norm[25] >= tolerance);
   0x000000000000024f <+591>:	movaps %xmm9,%xmm12
   0x0000000000000253 <+595>:	mulss  %xmm11,%xmm12
   0x0000000000000258 <+600>:	movss  %xmm12,0xf8(%rsp)

130	    norm_product[7][13] = (norm[7]*norm[29] >= tolerance);
   0x0000000000000262 <+610>:	movaps %xmm14,%xmm12
   0x0000000000000266 <+614>:	mulss  %xmm13,%xmm12
   0x000000000000026b <+619>:	movss  %xmm12,0xf0(%rsp)

131	    norm_product[4][2] = (norm[4]*norm[18] >= tolerance);
   0x0000000000000275 <+629>:	movaps %xmm6,%xmm12
   0x0000000000000279 <+633>:	mulss  %xmm10,%xmm12
   0x0000000000000283 <+643>:	movss  %xmm12,0x148(%rsp)

132	    norm_product[5][6] = (norm[5]*norm[22] >= tolerance);
   0x000000000000028d <+653>:	movaps %xmm8,%xmm12
   0x0000000000000291 <+657>:	mulss  %xmm1,%xmm12
   0x00000000000002ab <+683>:	movss  %xmm12,0xe8(%rsp)

133	    norm_product[6][10] = (norm[6]*norm[26] >= tolerance);
   0x00000000000002b5 <+693>:	movaps %xmm0,%xmm12
   0x00000000000002b9 <+697>:	mulss  %xmm11,%xmm12
   0x00000000000002d2 <+722>:	movss  %xmm12,0xe0(%rsp)

134	    norm_product[7][14] = (norm[7]*norm[30] >= tolerance);
   0x00000000000002dc <+732>:	movss  0x90(%rdx,%rsi,1),%xmm12
   0x00000000000002e6 <+742>:	mulss  %xmm13,%xmm12
   0x0000000000000301 <+769>:	movss  %xmm12,0xd8(%rsp)

135	    norm_product[4][3] = (norm[4]*norm[19] >= tolerance);
   0x000000000000027e <+638>:	mulss  %xmm5,%xmm10
   0x0000000000000296 <+662>:	movss  %xmm10,0x1d8(%rsp)

136	    norm_product[5][7] = (norm[5]*norm[23] >= tolerance);
   0x00000000000002a7 <+679>:	mulss  %xmm4,%xmm1
   0x00000000000002be <+702>:	movss  %xmm1,0xd0(%rsp)

137	    norm_product[6][11] = (norm[6]*norm[27] >= tolerance);
   0x00000000000002cd <+717>:	mulss  %xmm3,%xmm11
   0x00000000000002eb <+747>:	movss  %xmm11,0xc8(%rsp)

138	    norm_product[7][15] = (norm[7]*norm[31] >= tolerance);
   0x00000000000002fc <+764>:	mulss  %xmm2,%xmm13
   0x0000000000000317 <+791>:	movss  %xmm13,0xb8(%rsp)

139	    norm_product[8][0] = (norm[8]*norm[16] >= tolerance);
   0x00000000000002a0 <+672>:	movss  0x58(%rdx,%rsi,1),%xmm10
   0x00000000000002c7 <+711>:	movss  0x38(%rdx,%rsi,1),%xmm1
   0x0000000000000312 <+786>:	mulss  %xmm1,%xmm10
   0x000000000000032d <+813>:	movss  %xmm10,0x108(%rsp)

140	    norm_product[9][4] = (norm[9]*norm[20] >= tolerance);
   0x0000000000000321 <+801>:	movss  0x68(%rdx,%rsi,1),%xmm13
   0x0000000000000337 <+823>:	movss  0x3c(%rdx,%rsi,1),%xmm10
   0x000000000000033e <+830>:	mulss  %xmm10,%xmm13
   0x0000000000000354 <+852>:	movss  %xmm13,0xb0(%rsp)

141	    norm_product[10][8] = (norm[10]*norm[24] >= tolerance);
   0x00000000000002f5 <+757>:	movss  0x78(%rdx,%rsi,1),%xmm11
   0x000000000000030b <+779>:	movss  0x40(%rdx,%rsi,1),%xmm12
   0x0000000000000328 <+808>:	mulss  %xmm12,%xmm11
   0x0000000000000343 <+835>:	movss  %xmm11,0xa8(%rsp)

142	    norm_product[11][12] = (norm[11]*norm[28] >= tolerance);
   0x000000000000034d <+845>:	movss  0x44(%rdx,%rsi,1),%xmm11
   0x000000000000035e <+862>:	movss  0x88(%rdx,%rsi,1),%xmm13
   0x0000000000000368 <+872>:	mulss  %xmm11,%xmm13
   0x000000000000036d <+877>:	movss  %xmm13,0x98(%rsp)

143	    norm_product[8][1] = (norm[8]*norm[17] >= tolerance);
   0x0000000000000377 <+887>:	movaps %xmm7,%xmm13
   0x000000000000037b <+891>:	mulss  %xmm1,%xmm13
   0x0000000000000380 <+896>:	movss  %xmm13,0xc0(%rsp)

144	    norm_product[9][5] = (norm[9]*norm[21] >= tolerance);
   0x000000000000038a <+906>:	movaps %xmm15,%xmm13
   0x000000000000038e <+910>:	mulss  %xmm10,%xmm13
   0x0000000000000393 <+915>:	movss  %xmm13,0xa0(%rsp)

145	    norm_product[10][9] = (norm[10]*norm[25] >= tolerance);
   0x000000000000039d <+925>:	movaps %xmm9,%xmm13
   0x00000000000003a1 <+929>:	mulss  %xmm12,%xmm13
   0x00000000000003a6 <+934>:	movss  %xmm13,0x88(%rsp)

146	    norm_product[11][13] = (norm[11]*norm[29] >= tolerance);
   0x00000000000003b0 <+944>:	movaps %xmm14,%xmm13
   0x00000000000003b4 <+948>:	mulss  %xmm11,%xmm13
   0x00000000000003b9 <+953>:	movss  %xmm13,0x78(%rsp)

147	    norm_product[8][2] = (norm[8]*norm[18] >= tolerance);
   0x00000000000003c0 <+960>:	movaps %xmm6,%xmm13
   0x00000000000003c4 <+964>:	mulss  %xmm1,%xmm13
   0x00000000000003cd <+973>:	movss  %xmm13,0x128(%rsp)

148	    norm_product[9][6] = (norm[9]*norm[22] >= tolerance);
   0x00000000000003d7 <+983>:	movaps %xmm8,%xmm13
   0x00000000000003db <+987>:	mulss  %xmm10,%xmm13
   0x00000000000003f4 <+1012>:	movss  %xmm13,0x80(%rsp)

149	    norm_product[10][10] = (norm[10]*norm[26] >= tolerance);
   0x00000000000003fe <+1022>:	movaps %xmm0,%xmm13
   0x0000000000000402 <+1026>:	mulss  %xmm12,%xmm13
   0x000000000000041a <+1050>:	movss  %xmm13,0x60(%rsp)

150	    norm_product[11][14] = (norm[11]*norm[30] >= tolerance);
   0x0000000000000421 <+1057>:	movss  0x90(%rdx,%rsi,1),%xmm13
   0x000000000000042b <+1067>:	mulss  %xmm11,%xmm13
   0x0000000000000443 <+1091>:	movss  %xmm13,0x50(%rsp)

151	    norm_product[8][3] = (norm[8]*norm[19] >= tolerance);
   0x00000000000003c9 <+969>:	mulss  %xmm5,%xmm1
   0x00000000000003e0 <+992>:	movss  %xmm1,0x188(%rsp)

152	    norm_product[9][7] = (norm[9]*norm[23] >= tolerance);
   0x00000000000003ef <+1007>:	mulss  %xmm4,%xmm10
   0x0000000000000407 <+1031>:	movss  %xmm10,0x68(%rsp)

153	    norm_product[10][11] = (norm[10]*norm[27] >= tolerance);
   0x0000000000000415 <+1045>:	mulss  %xmm3,%xmm12
   0x0000000000000430 <+1072>:	movss  %xmm12,0x48(%rsp)

154	    norm_product[11][15] = (norm[11]*norm[31] >= tolerance);
   0x000000000000043e <+1086>:	mulss  %xmm2,%xmm11
   0x0000000000000456 <+1110>:	movss  %xmm11,0x38(%rsp)

155	    norm_product[12][0] = (norm[12]*norm[16] >= tolerance);
   0x00000000000003e9 <+1001>:	movss  0x48(%rdx,%rsi,1),%xmm1
   0x0000000000000437 <+1079>:	movss  0x58(%rdx,%rsi,1),%xmm12
   0x0000000000000451 <+1105>:	mulss  %xmm1,%xmm12
   0x0000000000000469 <+1129>:	movss  %xmm12,0x1e0(%rsp)

156	    norm_product[13][4] = (norm[13]*norm[20] >= tolerance);
   0x000000000000040e <+1038>:	movss  0x68(%rdx,%rsi,1),%xmm10
   0x000000000000044a <+1098>:	movss  0x4c(%rdx,%rsi,1),%xmm13
   0x0000000000000464 <+1124>:	mulss  %xmm13,%xmm10
   0x000000000000047f <+1151>:	movss  %xmm10,0x90(%rsp)

157	    norm_product[14][8] = (norm[14]*norm[24] >= tolerance);
   0x000000000000045d <+1117>:	movss  0x50(%rdx,%rsi,1),%xmm11
   0x0000000000000473 <+1139>:	movss  0x78(%rdx,%rsi,1),%xmm12
   0x000000000000047a <+1146>:	mulss  %xmm11,%xmm12
   0x0000000000000497 <+1175>:	movss  %xmm12,0x40(%rsp)

158	    norm_product[15][12] = (norm[15]*norm[28] >= tolerance);
   0x0000000000000489 <+1161>:	movss  0x88(%rdx,%rsi,1),%xmm10
   0x000000000000049e <+1182>:	movss  0x54(%rdx,%rsi,1),%xmm12
   0x00000000000004a5 <+1189>:	mulss  %xmm12,%xmm10
   0x00000000000004aa <+1194>:	movss  %xmm10,0x28(%rsp)

159	    norm_product[12][1] = (norm[12]*norm[17] >= tolerance);
   0x0000000000000493 <+1171>:	mulss  %xmm1,%xmm7

160	    norm_product[13][5] = (norm[13]*norm[21] >= tolerance);
   0x00000000000004bb <+1211>:	mulss  %xmm13,%xmm15
   0x00000000000004c0 <+1216>:	movss  %xmm15,0x70(%rsp)

161	    norm_product[14][9] = (norm[14]*norm[25] >= tolerance);
   0x00000000000004d1 <+1233>:	mulss  %xmm11,%xmm9
   0x00000000000004d6 <+1238>:	movss  %xmm9,0x30(%rsp)

162	    norm_product[15][13] = (norm[15]*norm[29] >= tolerance);
   0x00000000000004e7 <+1255>:	mulss  %xmm12,%xmm14
   0x00000000000004ec <+1260>:	movss  %xmm14,0x18(%rsp)

163	    norm_product[12][2] = (norm[12]*norm[18] >= tolerance);
   0x00000000000004f3 <+1267>:	mulss  %xmm1,%xmm6

164	    norm_product[13][6] = (norm[13]*norm[22] >= tolerance);
   0x00000000000004f7 <+1271>:	mulss  %xmm13,%xmm8
   0x00000000000004fc <+1276>:	movss  %xmm8,0x58(%rsp)

165	    norm_product[14][10] = (norm[14]*norm[26] >= tolerance);
   0x000000000000050d <+1293>:	mulss  %xmm11,%xmm0
   0x0000000000000512 <+1298>:	movss  %xmm0,0x20(%rsp)

166	    norm_product[15][14] = (norm[15]*norm[30] >= tolerance);
   0x00000000000004c7 <+1223>:	movss  0x90(%rdx,%rsi,1),%xmm15
   0x0000000000000521 <+1313>:	mulss  %xmm12,%xmm15
   0x0000000000000526 <+1318>:	movss  %xmm15,0x10(%rsp)

167	    norm_product[12][3] = (norm[12]*norm[19] >= tolerance);
   0x0000000000000531 <+1329>:	mulss  %xmm1,%xmm5

168	    norm_product[13][7] = (norm[13]*norm[23] >= tolerance);
   0x0000000000000545 <+1349>:	mulss  %xmm13,%xmm4

169	    norm_product[14][11] = (norm[14]*norm[27] >= tolerance);
   0x000000000000054e <+1358>:	mulss  %xmm11,%xmm3

170	    norm_product[15][15] = (norm[15]*norm[31] >= tolerance);
   0x0000000000000557 <+1367>:	mulss  %xmm12,%xmm2

171
172	    /* Reset C(1,1) matrix accumulators */
173	    C_row[0] = _mm_setzero_ps();
   0x000000000000055c <+1372>:	xorps  %xmm12,%xmm12

174	    C_row[1] = _mm_setzero_ps();
   0x000000000000052d <+1325>:	xorps  %xmm15,%xmm15

175	    C_row[2] = _mm_setzero_ps();
   0x0000000000000553 <+1363>:	xorps  %xmm11,%xmm11

176	    C_row[3] = _mm_setzero_ps();
   0x000000000000054a <+1354>:	xorps  %xmm13,%xmm13

177
178	    if (norm_product[0][0] &&
   0x00000000000004b1 <+1201>:	movss  0x1d0(%rsp),%xmm10
   0x00000000000004dd <+1245>:	movss  0x1d8(%rsp),%xmm9
   0x0000000000000503 <+1283>:	movss  0x1e0(%rsp),%xmm8
   0x0000000000000518 <+1304>:	movss  0x1c8(%rsp),%xmm0
   0x0000000000000535 <+1333>:	movss  0x1c0(%rsp),%xmm1
   0x000000000000053e <+1342>:	comiss %xmm1,%xmm0
   0x0000000000000541 <+1345>:	movaps (%rsp),%xmm0
   0x0000000000000560 <+1376>:	jb     0xa45 <stream_kernel+2629>

179	        norm_product[1][4] &&
   0x0000000000000566 <+1382>:	movss  0x1b0(%rsp),%xmm14
   0x0000000000000570 <+1392>:	comiss %xmm1,%xmm14
   0x0000000000000574 <+1396>:	jb     0xa45 <stream_kernel+2629>

180	        norm_product[2][8] &&
   0x000000000000057a <+1402>:	movss  0x1a8(%rsp),%xmm14
   0x0000000000000584 <+1412>:	comiss %xmm1,%xmm14
   0x0000000000000588 <+1416>:	jb     0xa45 <stream_kernel+2629>

181	        norm_product[3][12])
   0x000000000000058e <+1422>:	movss  0x1a0(%rsp),%xmm14
   0x0000000000000598 <+1432>:	comiss %xmm1,%xmm14
   0x000000000000059c <+1436>:	jb     0xa45 <stream_kernel+2629>

182	    {
183	      /* A(1,1)*B(1,1) = C(1,1). */
184	      for (i = 0; i < 4; i++)
185	      {
186	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
187	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
188	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000005a2 <+1442>:	movaps (%r8),%xmm12
   0x00000000000005ba <+1466>:	mulps  (%rcx),%xmm12
   0x00000000000005e4 <+1508>:	movaps 0x40(%r8),%xmm15
   0x00000000000005e9 <+1513>:	mulps  (%rcx),%xmm15
   0x00000000000005f1 <+1521>:	movaps 0x80(%r8),%xmm11
   0x00000000000005f9 <+1529>:	mulps  (%rcx),%xmm11
   0x0000000000000653 <+1619>:	movaps 0xc0(%r8),%xmm14
   0x000000000000065b <+1627>:	mulps  (%rcx),%xmm14

189
190	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
191	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
192	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000005a6 <+1446>:	movaps 0x10(%r8),%xmm15
   0x00000000000005b0 <+1456>:	movaps 0x50(%r8),%xmm13
   0x00000000000005be <+1470>:	mulps  0x10(%rcx),%xmm15
   0x00000000000005c8 <+1480>:	mulps  0x10(%rcx),%xmm13
   0x00000000000005d2 <+1490>:	addps  %xmm12,%xmm15
   0x00000000000005fd <+1533>:	addps  %xmm15,%xmm13
   0x000000000000060f <+1551>:	movaps 0x90(%r8),%xmm13
   0x0000000000000617 <+1559>:	mulps  0x10(%rcx),%xmm13
   0x000000000000062d <+1581>:	addps  %xmm11,%xmm13
   0x0000000000000642 <+1602>:	movaps 0xd0(%r8),%xmm13
   0x000000000000064a <+1610>:	mulps  0x10(%rcx),%xmm13
   0x000000000000065f <+1631>:	addps  %xmm14,%xmm13

193
194	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
195	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
196	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000005ab <+1451>:	movaps 0x20(%r8),%xmm11
   0x00000000000005b5 <+1461>:	movaps 0x60(%r8),%xmm14
   0x00000000000005c3 <+1475>:	mulps  0x20(%rcx),%xmm11
   0x00000000000005cd <+1485>:	mulps  0x20(%rcx),%xmm14
   0x00000000000005e0 <+1504>:	addps  %xmm15,%xmm11
   0x000000000000060b <+1547>:	addps  %xmm13,%xmm14
   0x0000000000000620 <+1568>:	movaps 0xa0(%r8),%xmm14
   0x0000000000000628 <+1576>:	mulps  0x20(%rcx),%xmm14
   0x000000000000063e <+1598>:	addps  %xmm13,%xmm14
   0x0000000000000663 <+1635>:	movaps 0xe0(%r8),%xmm14
   0x000000000000066b <+1643>:	mulps  0x20(%rcx),%xmm14
   0x0000000000000670 <+1648>:	addps  %xmm13,%xmm14

197
198	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
199	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
200	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000005d6 <+1494>:	movaps 0x30(%r8),%xmm12
   0x00000000000005db <+1499>:	mulps  0x30(%rcx),%xmm12
   0x00000000000005ed <+1517>:	addps  %xmm11,%xmm12
   0x0000000000000601 <+1537>:	movaps 0x70(%r8),%xmm15
   0x0000000000000606 <+1542>:	mulps  0x30(%rcx),%xmm15
   0x000000000000061c <+1564>:	addps  %xmm14,%xmm15
   0x0000000000000631 <+1585>:	movaps 0xb0(%r8),%xmm11
   0x0000000000000639 <+1593>:	mulps  0x30(%rcx),%xmm11
   0x000000000000064f <+1615>:	addps  %xmm14,%xmm11
   0x0000000000000674 <+1652>:	movaps 0xf0(%r8),%xmm13
   0x000000000000067c <+1660>:	mulps  0x30(%rcx),%xmm13
   0x0000000000000681 <+1665>:	addps  %xmm14,%xmm13

201	      }
202
203	      /* A(1,2)*B(2,1) = C(1,1). */
204	      for (i = 0; i < 4; i++)
205	      {
206	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
207	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
208	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000685 <+1669>:	movaps 0x100(%r8),%xmm14
   0x000000000000068d <+1677>:	mulps  0x100(%rcx),%xmm14
   0x0000000000000695 <+1685>:	addps  %xmm12,%xmm14
   0x00000000000006d5 <+1749>:	movaps 0x140(%r8),%xmm14
   0x00000000000006dd <+1757>:	mulps  0x100(%rcx),%xmm14
   0x00000000000006e5 <+1765>:	addps  %xmm15,%xmm14
   0x0000000000000725 <+1829>:	movaps 0x180(%r8),%xmm14
   0x000000000000072d <+1837>:	mulps  0x100(%rcx),%xmm14
   0x0000000000000735 <+1845>:	addps  %xmm11,%xmm14
   0x0000000000000775 <+1909>:	movaps 0x1c0(%r8),%xmm14
   0x000000000000077d <+1917>:	mulps  0x100(%rcx),%xmm14
   0x0000000000000785 <+1925>:	addps  %xmm13,%xmm14

209
210	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
211	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
212	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000699 <+1689>:	movaps 0x110(%r8),%xmm12
   0x00000000000006a1 <+1697>:	mulps  0x110(%rcx),%xmm12
   0x00000000000006a9 <+1705>:	addps  %xmm14,%xmm12
   0x00000000000006e9 <+1769>:	movaps 0x150(%r8),%xmm15
   0x00000000000006f1 <+1777>:	mulps  0x110(%rcx),%xmm15
   0x00000000000006f9 <+1785>:	addps  %xmm14,%xmm15
   0x0000000000000739 <+1849>:	movaps 0x190(%r8),%xmm11
   0x0000000000000741 <+1857>:	mulps  0x110(%rcx),%xmm11
   0x0000000000000749 <+1865>:	addps  %xmm14,%xmm11
   0x0000000000000789 <+1929>:	movaps 0x1d0(%r8),%xmm13
   0x0000000000000791 <+1937>:	mulps  0x110(%rcx),%xmm13
   0x0000000000000799 <+1945>:	addps  %xmm14,%xmm13

213
214	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
215	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
216	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000006ad <+1709>:	movaps 0x120(%r8),%xmm14
   0x00000000000006b5 <+1717>:	mulps  0x120(%rcx),%xmm14
   0x00000000000006bd <+1725>:	addps  %xmm12,%xmm14
   0x00000000000006fd <+1789>:	movaps 0x160(%r8),%xmm14
   0x0000000000000705 <+1797>:	mulps  0x120(%rcx),%xmm14
   0x000000000000070d <+1805>:	addps  %xmm15,%xmm14
   0x000000000000074d <+1869>:	movaps 0x1a0(%r8),%xmm14
   0x0000000000000755 <+1877>:	mulps  0x120(%rcx),%xmm14
   0x000000000000075d <+1885>:	addps  %xmm11,%xmm14
   0x000000000000079d <+1949>:	movaps 0x1e0(%r8),%xmm14
   0x00000000000007a5 <+1957>:	mulps  0x120(%rcx),%xmm14
   0x00000000000007ad <+1965>:	addps  %xmm13,%xmm14

217
218	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
219	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
220	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000006c1 <+1729>:	movaps 0x130(%r8),%xmm12
   0x00000000000006c9 <+1737>:	mulps  0x130(%rcx),%xmm12
   0x00000000000006d1 <+1745>:	addps  %xmm14,%xmm12
   0x0000000000000711 <+1809>:	movaps 0x170(%r8),%xmm15
   0x0000000000000719 <+1817>:	mulps  0x130(%rcx),%xmm15
   0x0000000000000721 <+1825>:	addps  %xmm14,%xmm15
   0x0000000000000761 <+1889>:	movaps 0x1b0(%r8),%xmm11
   0x0000000000000769 <+1897>:	mulps  0x130(%rcx),%xmm11
   0x0000000000000771 <+1905>:	addps  %xmm14,%xmm11
   0x00000000000007b1 <+1969>:	movaps 0x1f0(%r8),%xmm13
   0x00000000000007b9 <+1977>:	mulps  0x130(%rcx),%xmm13
   0x00000000000007c1 <+1985>:	addps  %xmm14,%xmm13

221	      }
222
223	      /* A(1,3)*B(3,1) = C(1,1). */
224	      for (i = 0; i < 4; i++)
225	      {
226	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
227	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
228	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000007c5 <+1989>:	movaps 0x200(%r8),%xmm14
   0x00000000000007cd <+1997>:	mulps  0x200(%rcx),%xmm14
   0x00000000000007d5 <+2005>:	addps  %xmm12,%xmm14
   0x0000000000000815 <+2069>:	movaps 0x240(%r8),%xmm14
   0x000000000000081d <+2077>:	mulps  0x200(%rcx),%xmm14
   0x0000000000000825 <+2085>:	addps  %xmm15,%xmm14
   0x0000000000000865 <+2149>:	movaps 0x280(%r8),%xmm14
   0x000000000000086d <+2157>:	mulps  0x200(%rcx),%xmm14
   0x0000000000000875 <+2165>:	addps  %xmm11,%xmm14
   0x00000000000008b5 <+2229>:	movaps 0x2c0(%r8),%xmm14
   0x00000000000008bd <+2237>:	mulps  0x200(%rcx),%xmm14
   0x00000000000008c5 <+2245>:	addps  %xmm13,%xmm14

229
230	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
231	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
232	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000007d9 <+2009>:	movaps 0x210(%r8),%xmm12
   0x00000000000007e1 <+2017>:	mulps  0x210(%rcx),%xmm12
   0x00000000000007e9 <+2025>:	addps  %xmm14,%xmm12
   0x0000000000000829 <+2089>:	movaps 0x250(%r8),%xmm15
   0x0000000000000831 <+2097>:	mulps  0x210(%rcx),%xmm15
   0x0000000000000839 <+2105>:	addps  %xmm14,%xmm15
   0x0000000000000879 <+2169>:	movaps 0x290(%r8),%xmm11
   0x0000000000000881 <+2177>:	mulps  0x210(%rcx),%xmm11
   0x0000000000000889 <+2185>:	addps  %xmm14,%xmm11
   0x00000000000008c9 <+2249>:	movaps 0x2d0(%r8),%xmm13
   0x00000000000008d1 <+2257>:	mulps  0x210(%rcx),%xmm13
   0x00000000000008d9 <+2265>:	addps  %xmm14,%xmm13

233
234	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
235	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
236	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000007ed <+2029>:	movaps 0x220(%r8),%xmm14
   0x00000000000007f5 <+2037>:	mulps  0x220(%rcx),%xmm14
   0x00000000000007fd <+2045>:	addps  %xmm12,%xmm14
   0x000000000000083d <+2109>:	movaps 0x260(%r8),%xmm14
   0x0000000000000845 <+2117>:	mulps  0x220(%rcx),%xmm14
   0x000000000000084d <+2125>:	addps  %xmm15,%xmm14
   0x000000000000088d <+2189>:	movaps 0x2a0(%r8),%xmm14
   0x0000000000000895 <+2197>:	mulps  0x220(%rcx),%xmm14
   0x000000000000089d <+2205>:	addps  %xmm11,%xmm14
   0x00000000000008dd <+2269>:	movaps 0x2e0(%r8),%xmm14
   0x00000000000008e5 <+2277>:	mulps  0x220(%rcx),%xmm14
   0x00000000000008ed <+2285>:	addps  %xmm13,%xmm14

237
238	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
239	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
240	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000801 <+2049>:	movaps 0x230(%r8),%xmm12
   0x0000000000000809 <+2057>:	mulps  0x230(%rcx),%xmm12
   0x0000000000000811 <+2065>:	addps  %xmm14,%xmm12
   0x0000000000000851 <+2129>:	movaps 0x270(%r8),%xmm15
   0x0000000000000859 <+2137>:	mulps  0x230(%rcx),%xmm15
   0x0000000000000861 <+2145>:	addps  %xmm14,%xmm15
   0x00000000000008a1 <+2209>:	movaps 0x2b0(%r8),%xmm11
   0x00000000000008a9 <+2217>:	mulps  0x230(%rcx),%xmm11
   0x00000000000008b1 <+2225>:	addps  %xmm14,%xmm11
   0x00000000000008f1 <+2289>:	movaps 0x2f0(%r8),%xmm13
   0x00000000000008f9 <+2297>:	mulps  0x230(%rcx),%xmm13
   0x0000000000000901 <+2305>:	addps  %xmm14,%xmm13

241	      }
242
243	      /* A(1,4)*B(4,1) = C(1,1). */
244	      for (i = 0; i < 4; i++)
245	      {
246	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
247	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
248	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000905 <+2309>:	movaps 0x300(%r8),%xmm14
   0x000000000000090d <+2317>:	mulps  0x300(%rcx),%xmm14
   0x0000000000000915 <+2325>:	addps  %xmm12,%xmm14
   0x0000000000000955 <+2389>:	movaps 0x340(%r8),%xmm14
   0x000000000000095d <+2397>:	mulps  0x300(%rcx),%xmm14
   0x0000000000000965 <+2405>:	addps  %xmm15,%xmm14
   0x00000000000009a5 <+2469>:	movaps 0x380(%r8),%xmm14
   0x00000000000009ad <+2477>:	mulps  0x300(%rcx),%xmm14
   0x00000000000009b5 <+2485>:	addps  %xmm11,%xmm14
   0x00000000000009f5 <+2549>:	movaps 0x3c0(%r8),%xmm14
   0x00000000000009fd <+2557>:	mulps  0x300(%rcx),%xmm14
   0x0000000000000a05 <+2565>:	addps  %xmm13,%xmm14

249
250	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
251	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
252	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000919 <+2329>:	movaps 0x310(%r8),%xmm12
   0x0000000000000921 <+2337>:	mulps  0x310(%rcx),%xmm12
   0x0000000000000929 <+2345>:	addps  %xmm14,%xmm12
   0x0000000000000969 <+2409>:	movaps 0x350(%r8),%xmm15
   0x0000000000000971 <+2417>:	mulps  0x310(%rcx),%xmm15
   0x0000000000000979 <+2425>:	addps  %xmm14,%xmm15
   0x00000000000009b9 <+2489>:	movaps 0x390(%r8),%xmm11
   0x00000000000009c1 <+2497>:	mulps  0x310(%rcx),%xmm11
   0x00000000000009c9 <+2505>:	addps  %xmm14,%xmm11
   0x0000000000000a09 <+2569>:	movaps 0x3d0(%r8),%xmm13
   0x0000000000000a11 <+2577>:	mulps  0x310(%rcx),%xmm13
   0x0000000000000a19 <+2585>:	addps  %xmm14,%xmm13

253
254	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
255	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
256	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000092d <+2349>:	movaps 0x320(%r8),%xmm14
   0x0000000000000935 <+2357>:	mulps  0x320(%rcx),%xmm14
   0x000000000000093d <+2365>:	addps  %xmm12,%xmm14
   0x000000000000097d <+2429>:	movaps 0x360(%r8),%xmm14
   0x0000000000000985 <+2437>:	mulps  0x320(%rcx),%xmm14
   0x000000000000098d <+2445>:	addps  %xmm15,%xmm14
   0x00000000000009cd <+2509>:	movaps 0x3a0(%r8),%xmm14
   0x00000000000009d5 <+2517>:	mulps  0x320(%rcx),%xmm14
   0x00000000000009dd <+2525>:	addps  %xmm11,%xmm14
   0x0000000000000a1d <+2589>:	movaps 0x3e0(%r8),%xmm14
   0x0000000000000a25 <+2597>:	mulps  0x320(%rcx),%xmm14
   0x0000000000000a2d <+2605>:	addps  %xmm13,%xmm14

257
258	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
259	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
260	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000941 <+2369>:	movaps 0x330(%r8),%xmm12
   0x0000000000000949 <+2377>:	mulps  0x330(%rcx),%xmm12
   0x0000000000000951 <+2385>:	addps  %xmm14,%xmm12
   0x0000000000000991 <+2449>:	movaps 0x370(%r8),%xmm15
   0x0000000000000999 <+2457>:	mulps  0x330(%rcx),%xmm15
   0x00000000000009a1 <+2465>:	addps  %xmm14,%xmm15
   0x00000000000009e1 <+2529>:	movaps 0x3b0(%r8),%xmm11
   0x00000000000009e9 <+2537>:	mulps  0x330(%rcx),%xmm11
   0x00000000000009f1 <+2545>:	addps  %xmm14,%xmm11
   0x0000000000000a31 <+2609>:	movaps 0x3f0(%r8),%xmm13
   0x0000000000000a39 <+2617>:	mulps  0x330(%rcx),%xmm13
   0x0000000000000a41 <+2625>:	addps  %xmm14,%xmm13

261	      }
262	    }
263
264	    /* Store C(1,1) block. */
265	    for (i = 0; i < 4; i++)
266	    {
267	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000000a45 <+2629>:	mulps  %xmm0,%xmm12
   0x0000000000000a5f <+2655>:	mulps  %xmm0,%xmm15
   0x0000000000000a71 <+2673>:	mulps  %xmm0,%xmm11
   0x0000000000000a83 <+2691>:	mulps  %xmm0,%xmm13

268	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_11]), C_row[i]);
   0x0000000000000a49 <+2633>:	addps  (%r9),%xmm12
   0x0000000000000a63 <+2659>:	addps  0x10(%r9),%xmm15
   0x0000000000000a75 <+2677>:	addps  0x20(%r9),%xmm11
   0x0000000000000a87 <+2695>:	addps  0x30(%r9),%xmm13

269	      _mm_store_ps(&C[i*4+C_OFFSET_11], C_row[i]);
   0x0000000000000a57 <+2647>:	movaps %xmm12,(%r9)
   0x0000000000000a68 <+2664>:	movaps %xmm15,0x10(%r9)
   0x0000000000000a7a <+2682>:	movaps %xmm11,0x20(%r9)
   0x0000000000000a8c <+2700>:	movaps %xmm13,0x30(%r9)

270	    }
271
272	    /* Reset C(1,2) matrix accumulators */
273	    C_row[0] = _mm_setzero_ps();
   0x0000000000000a5b <+2651>:	xorps  %xmm12,%xmm12

274	    C_row[1] = _mm_setzero_ps();
   0x0000000000000a6d <+2669>:	xorps  %xmm15,%xmm15

275	    C_row[2] = _mm_setzero_ps();
   0x0000000000000a7f <+2687>:	xorps  %xmm11,%xmm11

276	    C_row[3] = _mm_setzero_ps();
   0x0000000000000a91 <+2705>:	xorps  %xmm13,%xmm13

277
278	    if (norm_product[0][1] &&
   0x0000000000000a4d <+2637>:	movss  0x1b8(%rsp),%xmm14
   0x0000000000000a95 <+2709>:	comiss %xmm1,%xmm14
   0x0000000000000a99 <+2713>:	jb     0xf82 <stream_kernel+3970>

279	        norm_product[1][5] &&
   0x0000000000000a9f <+2719>:	movss  0x190(%rsp),%xmm14
   0x0000000000000aa9 <+2729>:	comiss %xmm1,%xmm14
   0x0000000000000aad <+2733>:	jb     0xf82 <stream_kernel+3970>

280	        norm_product[2][9] &&
   0x0000000000000ab3 <+2739>:	movss  0x180(%rsp),%xmm14
   0x0000000000000abd <+2749>:	comiss %xmm1,%xmm14
   0x0000000000000ac1 <+2753>:	jb     0xf82 <stream_kernel+3970>

281	        norm_product[3][13])
   0x0000000000000ac7 <+2759>:	movss  0x178(%rsp),%xmm14
   0x0000000000000ad1 <+2769>:	comiss %xmm1,%xmm14
   0x0000000000000ad5 <+2773>:	jb     0xf82 <stream_kernel+3970>

282	    {
283	      /* A(1,1)*B(1,2) = C(1,2). */
284	      for (i = 0; i < 4; i++)
285	      {
286	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
287	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
288	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000adb <+2779>:	movaps (%r8),%xmm12
   0x0000000000000af3 <+2803>:	mulps  0x40(%rcx),%xmm12
   0x0000000000000b1e <+2846>:	movaps 0x40(%r8),%xmm15
   0x0000000000000b23 <+2851>:	mulps  0x40(%rcx),%xmm15
   0x0000000000000b2c <+2860>:	movaps 0x80(%r8),%xmm11
   0x0000000000000b34 <+2868>:	mulps  0x40(%rcx),%xmm11
   0x0000000000000b8f <+2959>:	movaps 0xc0(%r8),%xmm14
   0x0000000000000b97 <+2967>:	mulps  0x40(%rcx),%xmm14

289
290	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
291	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
292	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000adf <+2783>:	movaps 0x10(%r8),%xmm15
   0x0000000000000ae9 <+2793>:	movaps 0x50(%r8),%xmm13
   0x0000000000000af8 <+2808>:	mulps  0x50(%rcx),%xmm15
   0x0000000000000b02 <+2818>:	mulps  0x50(%rcx),%xmm13
   0x0000000000000b0c <+2828>:	addps  %xmm12,%xmm15
   0x0000000000000b39 <+2873>:	addps  %xmm15,%xmm13
   0x0000000000000b4b <+2891>:	movaps 0x90(%r8),%xmm13
   0x0000000000000b53 <+2899>:	mulps  0x50(%rcx),%xmm13
   0x0000000000000b69 <+2921>:	addps  %xmm11,%xmm13
   0x0000000000000b7e <+2942>:	movaps 0xd0(%r8),%xmm13
   0x0000000000000b86 <+2950>:	mulps  0x50(%rcx),%xmm13
   0x0000000000000b9c <+2972>:	addps  %xmm14,%xmm13

293
294	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
295	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
296	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000ae4 <+2788>:	movaps 0x20(%r8),%xmm11
   0x0000000000000aee <+2798>:	movaps 0x60(%r8),%xmm14
   0x0000000000000afd <+2813>:	mulps  0x60(%rcx),%xmm11
   0x0000000000000b07 <+2823>:	mulps  0x60(%rcx),%xmm14
   0x0000000000000b1a <+2842>:	addps  %xmm15,%xmm11
   0x0000000000000b47 <+2887>:	addps  %xmm13,%xmm14
   0x0000000000000b5c <+2908>:	movaps 0xa0(%r8),%xmm14
   0x0000000000000b64 <+2916>:	mulps  0x60(%rcx),%xmm14
   0x0000000000000b7a <+2938>:	addps  %xmm13,%xmm14
   0x0000000000000ba0 <+2976>:	movaps 0xe0(%r8),%xmm14
   0x0000000000000ba8 <+2984>:	mulps  0x60(%rcx),%xmm14
   0x0000000000000bad <+2989>:	addps  %xmm13,%xmm14

297
298	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
299	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
300	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000b10 <+2832>:	movaps 0x30(%r8),%xmm12
   0x0000000000000b15 <+2837>:	mulps  0x70(%rcx),%xmm12
   0x0000000000000b28 <+2856>:	addps  %xmm11,%xmm12
   0x0000000000000b3d <+2877>:	movaps 0x70(%r8),%xmm15
   0x0000000000000b42 <+2882>:	mulps  0x70(%rcx),%xmm15
   0x0000000000000b58 <+2904>:	addps  %xmm14,%xmm15
   0x0000000000000b6d <+2925>:	movaps 0xb0(%r8),%xmm11
   0x0000000000000b75 <+2933>:	mulps  0x70(%rcx),%xmm11
   0x0000000000000b8b <+2955>:	addps  %xmm14,%xmm11
   0x0000000000000bb1 <+2993>:	movaps 0xf0(%r8),%xmm13
   0x0000000000000bb9 <+3001>:	mulps  0x70(%rcx),%xmm13
   0x0000000000000bbe <+3006>:	addps  %xmm14,%xmm13

301	      }
302
303	      /* A(1,2)*B(2,2) = C(1,2). */
304	      for (i = 0; i < 4; i++)
305	      {
306	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
307	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
308	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000bc2 <+3010>:	movaps 0x100(%r8),%xmm14
   0x0000000000000bca <+3018>:	mulps  0x140(%rcx),%xmm14
   0x0000000000000bd2 <+3026>:	addps  %xmm12,%xmm14
   0x0000000000000c12 <+3090>:	movaps 0x140(%r8),%xmm14
   0x0000000000000c1a <+3098>:	mulps  0x140(%rcx),%xmm14
   0x0000000000000c22 <+3106>:	addps  %xmm15,%xmm14
   0x0000000000000c62 <+3170>:	movaps 0x180(%r8),%xmm14
   0x0000000000000c6a <+3178>:	mulps  0x140(%rcx),%xmm14
   0x0000000000000c72 <+3186>:	addps  %xmm11,%xmm14
   0x0000000000000cb2 <+3250>:	movaps 0x1c0(%r8),%xmm14
   0x0000000000000cba <+3258>:	mulps  0x140(%rcx),%xmm14
   0x0000000000000cc2 <+3266>:	addps  %xmm13,%xmm14

309
310	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
311	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
312	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000bd6 <+3030>:	movaps 0x110(%r8),%xmm12
   0x0000000000000bde <+3038>:	mulps  0x150(%rcx),%xmm12
   0x0000000000000be6 <+3046>:	addps  %xmm14,%xmm12
   0x0000000000000c26 <+3110>:	movaps 0x150(%r8),%xmm15
   0x0000000000000c2e <+3118>:	mulps  0x150(%rcx),%xmm15
   0x0000000000000c36 <+3126>:	addps  %xmm14,%xmm15
   0x0000000000000c76 <+3190>:	movaps 0x190(%r8),%xmm11
   0x0000000000000c7e <+3198>:	mulps  0x150(%rcx),%xmm11
   0x0000000000000c86 <+3206>:	addps  %xmm14,%xmm11
   0x0000000000000cc6 <+3270>:	movaps 0x1d0(%r8),%xmm13
   0x0000000000000cce <+3278>:	mulps  0x150(%rcx),%xmm13
   0x0000000000000cd6 <+3286>:	addps  %xmm14,%xmm13

313
314	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
315	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
316	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000bea <+3050>:	movaps 0x120(%r8),%xmm14
   0x0000000000000bf2 <+3058>:	mulps  0x160(%rcx),%xmm14
   0x0000000000000bfa <+3066>:	addps  %xmm12,%xmm14
   0x0000000000000c3a <+3130>:	movaps 0x160(%r8),%xmm14
   0x0000000000000c42 <+3138>:	mulps  0x160(%rcx),%xmm14
   0x0000000000000c4a <+3146>:	addps  %xmm15,%xmm14
   0x0000000000000c8a <+3210>:	movaps 0x1a0(%r8),%xmm14
   0x0000000000000c92 <+3218>:	mulps  0x160(%rcx),%xmm14
   0x0000000000000c9a <+3226>:	addps  %xmm11,%xmm14
   0x0000000000000cda <+3290>:	movaps 0x1e0(%r8),%xmm14
   0x0000000000000ce2 <+3298>:	mulps  0x160(%rcx),%xmm14
   0x0000000000000cea <+3306>:	addps  %xmm13,%xmm14

317
318	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
319	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
320	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000bfe <+3070>:	movaps 0x130(%r8),%xmm12
   0x0000000000000c06 <+3078>:	mulps  0x170(%rcx),%xmm12
   0x0000000000000c0e <+3086>:	addps  %xmm14,%xmm12
   0x0000000000000c4e <+3150>:	movaps 0x170(%r8),%xmm15
   0x0000000000000c56 <+3158>:	mulps  0x170(%rcx),%xmm15
   0x0000000000000c5e <+3166>:	addps  %xmm14,%xmm15
   0x0000000000000c9e <+3230>:	movaps 0x1b0(%r8),%xmm11
   0x0000000000000ca6 <+3238>:	mulps  0x170(%rcx),%xmm11
   0x0000000000000cae <+3246>:	addps  %xmm14,%xmm11
   0x0000000000000cee <+3310>:	movaps 0x1f0(%r8),%xmm13
   0x0000000000000cf6 <+3318>:	mulps  0x170(%rcx),%xmm13
   0x0000000000000cfe <+3326>:	addps  %xmm14,%xmm13

321	      }
322
323	      /* A(1,3)*B(3,2) = C(1,2). */
324	      for (i = 0; i < 4; i++)
325	      {
326	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
327	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
328	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d02 <+3330>:	movaps 0x200(%r8),%xmm14
   0x0000000000000d0a <+3338>:	mulps  0x240(%rcx),%xmm14
   0x0000000000000d12 <+3346>:	addps  %xmm12,%xmm14
   0x0000000000000d52 <+3410>:	movaps 0x240(%r8),%xmm14
   0x0000000000000d5a <+3418>:	mulps  0x240(%rcx),%xmm14
   0x0000000000000d62 <+3426>:	addps  %xmm15,%xmm14
   0x0000000000000da2 <+3490>:	movaps 0x280(%r8),%xmm14
   0x0000000000000daa <+3498>:	mulps  0x240(%rcx),%xmm14
   0x0000000000000db2 <+3506>:	addps  %xmm11,%xmm14
   0x0000000000000df2 <+3570>:	movaps 0x2c0(%r8),%xmm14
   0x0000000000000dfa <+3578>:	mulps  0x240(%rcx),%xmm14
   0x0000000000000e02 <+3586>:	addps  %xmm13,%xmm14

329
330	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
331	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
332	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d16 <+3350>:	movaps 0x210(%r8),%xmm12
   0x0000000000000d1e <+3358>:	mulps  0x250(%rcx),%xmm12
   0x0000000000000d26 <+3366>:	addps  %xmm14,%xmm12
   0x0000000000000d66 <+3430>:	movaps 0x250(%r8),%xmm15
   0x0000000000000d6e <+3438>:	mulps  0x250(%rcx),%xmm15
   0x0000000000000d76 <+3446>:	addps  %xmm14,%xmm15
   0x0000000000000db6 <+3510>:	movaps 0x290(%r8),%xmm11
   0x0000000000000dbe <+3518>:	mulps  0x250(%rcx),%xmm11
   0x0000000000000dc6 <+3526>:	addps  %xmm14,%xmm11
   0x0000000000000e06 <+3590>:	movaps 0x2d0(%r8),%xmm13
   0x0000000000000e0e <+3598>:	mulps  0x250(%rcx),%xmm13
   0x0000000000000e16 <+3606>:	addps  %xmm14,%xmm13

333
334	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
335	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
336	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d2a <+3370>:	movaps 0x220(%r8),%xmm14
   0x0000000000000d32 <+3378>:	mulps  0x260(%rcx),%xmm14
   0x0000000000000d3a <+3386>:	addps  %xmm12,%xmm14
   0x0000000000000d7a <+3450>:	movaps 0x260(%r8),%xmm14
   0x0000000000000d82 <+3458>:	mulps  0x260(%rcx),%xmm14
   0x0000000000000d8a <+3466>:	addps  %xmm15,%xmm14
   0x0000000000000dca <+3530>:	movaps 0x2a0(%r8),%xmm14
   0x0000000000000dd2 <+3538>:	mulps  0x260(%rcx),%xmm14
   0x0000000000000dda <+3546>:	addps  %xmm11,%xmm14
   0x0000000000000e1a <+3610>:	movaps 0x2e0(%r8),%xmm14
   0x0000000000000e22 <+3618>:	mulps  0x260(%rcx),%xmm14
   0x0000000000000e2a <+3626>:	addps  %xmm13,%xmm14

337
338	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
339	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
340	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d3e <+3390>:	movaps 0x230(%r8),%xmm12
   0x0000000000000d46 <+3398>:	mulps  0x270(%rcx),%xmm12
   0x0000000000000d4e <+3406>:	addps  %xmm14,%xmm12
   0x0000000000000d8e <+3470>:	movaps 0x270(%r8),%xmm15
   0x0000000000000d96 <+3478>:	mulps  0x270(%rcx),%xmm15
   0x0000000000000d9e <+3486>:	addps  %xmm14,%xmm15
   0x0000000000000dde <+3550>:	movaps 0x2b0(%r8),%xmm11
   0x0000000000000de6 <+3558>:	mulps  0x270(%rcx),%xmm11
   0x0000000000000dee <+3566>:	addps  %xmm14,%xmm11
   0x0000000000000e2e <+3630>:	movaps 0x2f0(%r8),%xmm13
   0x0000000000000e36 <+3638>:	mulps  0x270(%rcx),%xmm13
   0x0000000000000e3e <+3646>:	addps  %xmm14,%xmm13

341	      }
342
343	      /* A(1,4)*B(4,2) = C(1,2). */
344	      for (i = 0; i < 4; i++)
345	      {
346	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
347	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
348	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000e42 <+3650>:	movaps 0x300(%r8),%xmm14
   0x0000000000000e4a <+3658>:	mulps  0x340(%rcx),%xmm14
   0x0000000000000e52 <+3666>:	addps  %xmm12,%xmm14
   0x0000000000000e92 <+3730>:	movaps 0x340(%r8),%xmm14
   0x0000000000000e9a <+3738>:	mulps  0x340(%rcx),%xmm14
   0x0000000000000ea2 <+3746>:	addps  %xmm15,%xmm14
   0x0000000000000ee2 <+3810>:	movaps 0x380(%r8),%xmm14
   0x0000000000000eea <+3818>:	mulps  0x340(%rcx),%xmm14
   0x0000000000000ef2 <+3826>:	addps  %xmm11,%xmm14
   0x0000000000000f32 <+3890>:	movaps 0x3c0(%r8),%xmm14
   0x0000000000000f3a <+3898>:	mulps  0x340(%rcx),%xmm14
   0x0000000000000f42 <+3906>:	addps  %xmm13,%xmm14

349
350	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
351	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
352	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000e56 <+3670>:	movaps 0x310(%r8),%xmm12
   0x0000000000000e5e <+3678>:	mulps  0x350(%rcx),%xmm12
   0x0000000000000e66 <+3686>:	addps  %xmm14,%xmm12
   0x0000000000000ea6 <+3750>:	movaps 0x350(%r8),%xmm15
   0x0000000000000eae <+3758>:	mulps  0x350(%rcx),%xmm15
   0x0000000000000eb6 <+3766>:	addps  %xmm14,%xmm15
   0x0000000000000ef6 <+3830>:	movaps 0x390(%r8),%xmm11
   0x0000000000000efe <+3838>:	mulps  0x350(%rcx),%xmm11
   0x0000000000000f06 <+3846>:	addps  %xmm14,%xmm11
   0x0000000000000f46 <+3910>:	movaps 0x3d0(%r8),%xmm13
   0x0000000000000f4e <+3918>:	mulps  0x350(%rcx),%xmm13
   0x0000000000000f56 <+3926>:	addps  %xmm14,%xmm13

353
354	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
355	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
356	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000e6a <+3690>:	movaps 0x320(%r8),%xmm14
   0x0000000000000e72 <+3698>:	mulps  0x360(%rcx),%xmm14
   0x0000000000000e7a <+3706>:	addps  %xmm12,%xmm14
   0x0000000000000eba <+3770>:	movaps 0x360(%r8),%xmm14
   0x0000000000000ec2 <+3778>:	mulps  0x360(%rcx),%xmm14
   0x0000000000000eca <+3786>:	addps  %xmm15,%xmm14
   0x0000000000000f0a <+3850>:	movaps 0x3a0(%r8),%xmm14
   0x0000000000000f12 <+3858>:	mulps  0x360(%rcx),%xmm14
   0x0000000000000f1a <+3866>:	addps  %xmm11,%xmm14
   0x0000000000000f5a <+3930>:	movaps 0x3e0(%r8),%xmm14
   0x0000000000000f62 <+3938>:	mulps  0x360(%rcx),%xmm14
   0x0000000000000f6a <+3946>:	addps  %xmm13,%xmm14

357
358	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
359	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
360	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000e7e <+3710>:	movaps 0x330(%r8),%xmm12
   0x0000000000000e86 <+3718>:	mulps  0x370(%rcx),%xmm12
   0x0000000000000e8e <+3726>:	addps  %xmm14,%xmm12
   0x0000000000000ece <+3790>:	movaps 0x370(%r8),%xmm15
   0x0000000000000ed6 <+3798>:	mulps  0x370(%rcx),%xmm15
   0x0000000000000ede <+3806>:	addps  %xmm14,%xmm15
   0x0000000000000f1e <+3870>:	movaps 0x3b0(%r8),%xmm11
   0x0000000000000f26 <+3878>:	mulps  0x370(%rcx),%xmm11
   0x0000000000000f2e <+3886>:	addps  %xmm14,%xmm11
   0x0000000000000f6e <+3950>:	movaps 0x3f0(%r8),%xmm13
   0x0000000000000f76 <+3958>:	mulps  0x370(%rcx),%xmm13
   0x0000000000000f7e <+3966>:	addps  %xmm14,%xmm13

361	      }
362	    }
363
364	    /* Store C(1,2) block. */
365	    for (i = 0; i < 4; i++)
366	    {
367	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000000f82 <+3970>:	mulps  %xmm0,%xmm12
   0x0000000000000f9e <+3998>:	mulps  %xmm0,%xmm15
   0x0000000000000fb0 <+4016>:	mulps  %xmm0,%xmm11
   0x0000000000000fc2 <+4034>:	mulps  %xmm0,%xmm13

368	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_12]), C_row[i]);
   0x0000000000000f86 <+3974>:	addps  0x40(%r9),%xmm12
   0x0000000000000fa2 <+4002>:	addps  0x50(%r9),%xmm15
   0x0000000000000fb4 <+4020>:	addps  0x60(%r9),%xmm11
   0x0000000000000fc6 <+4038>:	addps  0x70(%r9),%xmm13

369	      _mm_store_ps(&C[i*4+C_OFFSET_12], C_row[i]);
   0x0000000000000f95 <+3989>:	movaps %xmm12,0x40(%r9)
   0x0000000000000fa7 <+4007>:	movaps %xmm15,0x50(%r9)
   0x0000000000000fb9 <+4025>:	movaps %xmm11,0x60(%r9)
   0x0000000000000fcb <+4043>:	movaps %xmm13,0x70(%r9)

370	    }
371
372	    /* Reset C(1,3) matrix accumulators */
373	    C_row[0] = _mm_setzero_ps();
   0x0000000000000f9a <+3994>:	xorps  %xmm12,%xmm12

374	    C_row[1] = _mm_setzero_ps();
   0x0000000000000fac <+4012>:	xorps  %xmm15,%xmm15

375	    C_row[2] = _mm_setzero_ps();
   0x0000000000000fbe <+4030>:	xorps  %xmm11,%xmm11

376	    C_row[3] = _mm_setzero_ps();
   0x0000000000000fd0 <+4048>:	xorps  %xmm13,%xmm13

377
378	    if (norm_product[0][2] &&
   0x0000000000000f8b <+3979>:	movss  0x198(%rsp),%xmm14
   0x0000000000000fd4 <+4052>:	comiss %xmm1,%xmm14
   0x0000000000000fd8 <+4056>:	jb     0x14f1 <stream_kernel+5361>

379	        norm_product[1][6] &&
   0x0000000000000fde <+4062>:	movss  0x168(%rsp),%xmm14
   0x0000000000000fe8 <+4072>:	comiss %xmm1,%xmm14
   0x0000000000000fec <+4076>:	jb     0x14f1 <stream_kernel+5361>

380	        norm_product[2][10] &&
   0x0000000000000ff2 <+4082>:	movss  0x158(%rsp),%xmm14
   0x0000000000000ffc <+4092>:	comiss %xmm1,%xmm14
   0x0000000000001000 <+4096>:	jb     0x14f1 <stream_kernel+5361>

381	        norm_product[3][14])
   0x0000000000001006 <+4102>:	movss  0x150(%rsp),%xmm14
   0x0000000000001010 <+4112>:	comiss %xmm1,%xmm14
   0x0000000000001014 <+4116>:	jb     0x14f1 <stream_kernel+5361>

382	    {
383	      /* A(1,1)*B(1,3) = C(1,3). */
384	      for (i = 0; i < 4; i++)
385	      {
386	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
387	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
388	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000101a <+4122>:	movaps (%r8),%xmm12
   0x0000000000001032 <+4146>:	mulps  0x80(%rcx),%xmm12
   0x000000000000106f <+4207>:	movaps 0x40(%r8),%xmm15
   0x0000000000001074 <+4212>:	mulps  0x80(%rcx),%xmm15
   0x0000000000001080 <+4224>:	movaps 0x80(%r8),%xmm11
   0x0000000000001088 <+4232>:	mulps  0x80(%rcx),%xmm11
   0x00000000000010f5 <+4341>:	movaps 0xc0(%r8),%xmm14
   0x00000000000010fd <+4349>:	mulps  0x80(%rcx),%xmm14

389
390	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
391	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
392	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000101e <+4126>:	movaps 0x10(%r8),%xmm15
   0x0000000000001028 <+4136>:	movaps 0x50(%r8),%xmm13
   0x000000000000103a <+4154>:	mulps  0x90(%rcx),%xmm15
   0x000000000000104a <+4170>:	mulps  0x90(%rcx),%xmm13
   0x000000000000105a <+4186>:	addps  %xmm12,%xmm15
   0x0000000000001090 <+4240>:	addps  %xmm15,%xmm13
   0x00000000000010a5 <+4261>:	movaps 0x90(%r8),%xmm13
   0x00000000000010ad <+4269>:	mulps  0x90(%rcx),%xmm13
   0x00000000000010c9 <+4297>:	addps  %xmm11,%xmm13
   0x00000000000010e1 <+4321>:	movaps 0xd0(%r8),%xmm13
   0x00000000000010e9 <+4329>:	mulps  0x90(%rcx),%xmm13
   0x0000000000001105 <+4357>:	addps  %xmm14,%xmm13

393
394	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
395	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
396	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001023 <+4131>:	movaps 0x20(%r8),%xmm11
   0x000000000000102d <+4141>:	movaps 0x60(%r8),%xmm14
   0x0000000000001042 <+4162>:	mulps  0xa0(%rcx),%xmm11
   0x0000000000001052 <+4178>:	mulps  0xa0(%rcx),%xmm14
   0x000000000000106b <+4203>:	addps  %xmm15,%xmm11
   0x00000000000010a1 <+4257>:	addps  %xmm13,%xmm14
   0x00000000000010b9 <+4281>:	movaps 0xa0(%r8),%xmm14
   0x00000000000010c1 <+4289>:	mulps  0xa0(%rcx),%xmm14
   0x00000000000010dd <+4317>:	addps  %xmm13,%xmm14
   0x0000000000001109 <+4361>:	movaps 0xe0(%r8),%xmm14
   0x0000000000001111 <+4369>:	mulps  0xa0(%rcx),%xmm14
   0x0000000000001119 <+4377>:	addps  %xmm13,%xmm14

397
398	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
399	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
400	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000105e <+4190>:	movaps 0x30(%r8),%xmm12
   0x0000000000001063 <+4195>:	mulps  0xb0(%rcx),%xmm12
   0x000000000000107c <+4220>:	addps  %xmm11,%xmm12
   0x0000000000001094 <+4244>:	movaps 0x70(%r8),%xmm15
   0x0000000000001099 <+4249>:	mulps  0xb0(%rcx),%xmm15
   0x00000000000010b5 <+4277>:	addps  %xmm14,%xmm15
   0x00000000000010cd <+4301>:	movaps 0xb0(%r8),%xmm11
   0x00000000000010d5 <+4309>:	mulps  0xb0(%rcx),%xmm11
   0x00000000000010f1 <+4337>:	addps  %xmm14,%xmm11
   0x000000000000111d <+4381>:	movaps 0xf0(%r8),%xmm13
   0x0000000000001125 <+4389>:	mulps  0xb0(%rcx),%xmm13
   0x000000000000112d <+4397>:	addps  %xmm14,%xmm13

401	      }
402
403	      /* A(1,2)*B(2,3) = C(1,3). */
404	      for (i = 0; i < 4; i++)
405	      {
406	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
407	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
408	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001131 <+4401>:	movaps 0x100(%r8),%xmm14
   0x0000000000001139 <+4409>:	mulps  0x180(%rcx),%xmm14
   0x0000000000001141 <+4417>:	addps  %xmm12,%xmm14
   0x0000000000001181 <+4481>:	movaps 0x140(%r8),%xmm14
   0x0000000000001189 <+4489>:	mulps  0x180(%rcx),%xmm14
   0x0000000000001191 <+4497>:	addps  %xmm15,%xmm14
   0x00000000000011d1 <+4561>:	movaps 0x180(%r8),%xmm14
   0x00000000000011d9 <+4569>:	mulps  0x180(%rcx),%xmm14
   0x00000000000011e1 <+4577>:	addps  %xmm11,%xmm14
   0x0000000000001221 <+4641>:	movaps 0x1c0(%r8),%xmm14
   0x0000000000001229 <+4649>:	mulps  0x180(%rcx),%xmm14
   0x0000000000001231 <+4657>:	addps  %xmm13,%xmm14

409
410	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
411	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
412	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001145 <+4421>:	movaps 0x110(%r8),%xmm12
   0x000000000000114d <+4429>:	mulps  0x190(%rcx),%xmm12
   0x0000000000001155 <+4437>:	addps  %xmm14,%xmm12
   0x0000000000001195 <+4501>:	movaps 0x150(%r8),%xmm15
   0x000000000000119d <+4509>:	mulps  0x190(%rcx),%xmm15
   0x00000000000011a5 <+4517>:	addps  %xmm14,%xmm15
   0x00000000000011e5 <+4581>:	movaps 0x190(%r8),%xmm11
   0x00000000000011ed <+4589>:	mulps  0x190(%rcx),%xmm11
   0x00000000000011f5 <+4597>:	addps  %xmm14,%xmm11
   0x0000000000001235 <+4661>:	movaps 0x1d0(%r8),%xmm13
   0x000000000000123d <+4669>:	mulps  0x190(%rcx),%xmm13
   0x0000000000001245 <+4677>:	addps  %xmm14,%xmm13

413
414	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
415	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
416	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001159 <+4441>:	movaps 0x120(%r8),%xmm14
   0x0000000000001161 <+4449>:	mulps  0x1a0(%rcx),%xmm14
   0x0000000000001169 <+4457>:	addps  %xmm12,%xmm14
   0x00000000000011a9 <+4521>:	movaps 0x160(%r8),%xmm14
   0x00000000000011b1 <+4529>:	mulps  0x1a0(%rcx),%xmm14
   0x00000000000011b9 <+4537>:	addps  %xmm15,%xmm14
   0x00000000000011f9 <+4601>:	movaps 0x1a0(%r8),%xmm14
   0x0000000000001201 <+4609>:	mulps  0x1a0(%rcx),%xmm14
   0x0000000000001209 <+4617>:	addps  %xmm11,%xmm14
   0x0000000000001249 <+4681>:	movaps 0x1e0(%r8),%xmm14
   0x0000000000001251 <+4689>:	mulps  0x1a0(%rcx),%xmm14
   0x0000000000001259 <+4697>:	addps  %xmm13,%xmm14

417
418	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
419	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
420	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000116d <+4461>:	movaps 0x130(%r8),%xmm12
   0x0000000000001175 <+4469>:	mulps  0x1b0(%rcx),%xmm12
   0x000000000000117d <+4477>:	addps  %xmm14,%xmm12
   0x00000000000011bd <+4541>:	movaps 0x170(%r8),%xmm15
   0x00000000000011c5 <+4549>:	mulps  0x1b0(%rcx),%xmm15
   0x00000000000011cd <+4557>:	addps  %xmm14,%xmm15
   0x000000000000120d <+4621>:	movaps 0x1b0(%r8),%xmm11
   0x0000000000001215 <+4629>:	mulps  0x1b0(%rcx),%xmm11
   0x000000000000121d <+4637>:	addps  %xmm14,%xmm11
   0x000000000000125d <+4701>:	movaps 0x1f0(%r8),%xmm13
   0x0000000000001265 <+4709>:	mulps  0x1b0(%rcx),%xmm13
   0x000000000000126d <+4717>:	addps  %xmm14,%xmm13

421	      }
422
423	      /* A(1,3)*B(3,3) = C(1,3). */
424	      for (i = 0; i < 4; i++)
425	      {
426	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
427	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
428	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001271 <+4721>:	movaps 0x200(%r8),%xmm14
   0x0000000000001279 <+4729>:	mulps  0x280(%rcx),%xmm14
   0x0000000000001281 <+4737>:	addps  %xmm12,%xmm14
   0x00000000000012c1 <+4801>:	movaps 0x240(%r8),%xmm14
   0x00000000000012c9 <+4809>:	mulps  0x280(%rcx),%xmm14
   0x00000000000012d1 <+4817>:	addps  %xmm15,%xmm14
   0x0000000000001311 <+4881>:	movaps 0x280(%r8),%xmm14
   0x0000000000001319 <+4889>:	mulps  0x280(%rcx),%xmm14
   0x0000000000001321 <+4897>:	addps  %xmm11,%xmm14
   0x0000000000001361 <+4961>:	movaps 0x2c0(%r8),%xmm14
   0x0000000000001369 <+4969>:	mulps  0x280(%rcx),%xmm14
   0x0000000000001371 <+4977>:	addps  %xmm13,%xmm14

429
430	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
431	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
432	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001285 <+4741>:	movaps 0x210(%r8),%xmm12
   0x000000000000128d <+4749>:	mulps  0x290(%rcx),%xmm12
   0x0000000000001295 <+4757>:	addps  %xmm14,%xmm12
   0x00000000000012d5 <+4821>:	movaps 0x250(%r8),%xmm15
   0x00000000000012dd <+4829>:	mulps  0x290(%rcx),%xmm15
   0x00000000000012e5 <+4837>:	addps  %xmm14,%xmm15
   0x0000000000001325 <+4901>:	movaps 0x290(%r8),%xmm11
   0x000000000000132d <+4909>:	mulps  0x290(%rcx),%xmm11
   0x0000000000001335 <+4917>:	addps  %xmm14,%xmm11
   0x0000000000001375 <+4981>:	movaps 0x2d0(%r8),%xmm13
   0x000000000000137d <+4989>:	mulps  0x290(%rcx),%xmm13
   0x0000000000001385 <+4997>:	addps  %xmm14,%xmm13

433
434	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
435	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
436	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001299 <+4761>:	movaps 0x220(%r8),%xmm14
   0x00000000000012a1 <+4769>:	mulps  0x2a0(%rcx),%xmm14
   0x00000000000012a9 <+4777>:	addps  %xmm12,%xmm14
   0x00000000000012e9 <+4841>:	movaps 0x260(%r8),%xmm14
   0x00000000000012f1 <+4849>:	mulps  0x2a0(%rcx),%xmm14
   0x00000000000012f9 <+4857>:	addps  %xmm15,%xmm14
   0x0000000000001339 <+4921>:	movaps 0x2a0(%r8),%xmm14
   0x0000000000001341 <+4929>:	mulps  0x2a0(%rcx),%xmm14
   0x0000000000001349 <+4937>:	addps  %xmm11,%xmm14
   0x0000000000001389 <+5001>:	movaps 0x2e0(%r8),%xmm14
   0x0000000000001391 <+5009>:	mulps  0x2a0(%rcx),%xmm14
   0x0000000000001399 <+5017>:	addps  %xmm13,%xmm14

437
438	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
439	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
440	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000012ad <+4781>:	movaps 0x230(%r8),%xmm12
   0x00000000000012b5 <+4789>:	mulps  0x2b0(%rcx),%xmm12
   0x00000000000012bd <+4797>:	addps  %xmm14,%xmm12
   0x00000000000012fd <+4861>:	movaps 0x270(%r8),%xmm15
   0x0000000000001305 <+4869>:	mulps  0x2b0(%rcx),%xmm15
   0x000000000000130d <+4877>:	addps  %xmm14,%xmm15
   0x000000000000134d <+4941>:	movaps 0x2b0(%r8),%xmm11
   0x0000000000001355 <+4949>:	mulps  0x2b0(%rcx),%xmm11
   0x000000000000135d <+4957>:	addps  %xmm14,%xmm11
   0x000000000000139d <+5021>:	movaps 0x2f0(%r8),%xmm13
   0x00000000000013a5 <+5029>:	mulps  0x2b0(%rcx),%xmm13
   0x00000000000013ad <+5037>:	addps  %xmm14,%xmm13

441	      }
442
443	      /* A(1,4)*B(4,3) = C(1,3). */
444	      for (i = 0; i < 4; i++)
445	      {
446	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
447	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
448	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000013b1 <+5041>:	movaps 0x300(%r8),%xmm14
   0x00000000000013b9 <+5049>:	mulps  0x380(%rcx),%xmm14
   0x00000000000013c1 <+5057>:	addps  %xmm12,%xmm14
   0x0000000000001401 <+5121>:	movaps 0x340(%r8),%xmm14
   0x0000000000001409 <+5129>:	mulps  0x380(%rcx),%xmm14
   0x0000000000001411 <+5137>:	addps  %xmm15,%xmm14
   0x0000000000001451 <+5201>:	movaps 0x380(%r8),%xmm14
   0x0000000000001459 <+5209>:	mulps  0x380(%rcx),%xmm14
   0x0000000000001461 <+5217>:	addps  %xmm11,%xmm14
   0x00000000000014a1 <+5281>:	movaps 0x3c0(%r8),%xmm14
   0x00000000000014a9 <+5289>:	mulps  0x380(%rcx),%xmm14
   0x00000000000014b1 <+5297>:	addps  %xmm13,%xmm14

449
450	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
451	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
452	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000013c5 <+5061>:	movaps 0x310(%r8),%xmm12
   0x00000000000013cd <+5069>:	mulps  0x390(%rcx),%xmm12
   0x00000000000013d5 <+5077>:	addps  %xmm14,%xmm12
   0x0000000000001415 <+5141>:	movaps 0x350(%r8),%xmm15
   0x000000000000141d <+5149>:	mulps  0x390(%rcx),%xmm15
   0x0000000000001425 <+5157>:	addps  %xmm14,%xmm15
   0x0000000000001465 <+5221>:	movaps 0x390(%r8),%xmm11
   0x000000000000146d <+5229>:	mulps  0x390(%rcx),%xmm11
   0x0000000000001475 <+5237>:	addps  %xmm14,%xmm11
   0x00000000000014b5 <+5301>:	movaps 0x3d0(%r8),%xmm13
   0x00000000000014bd <+5309>:	mulps  0x390(%rcx),%xmm13
   0x00000000000014c5 <+5317>:	addps  %xmm14,%xmm13

453
454	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
455	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
456	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000013d9 <+5081>:	movaps 0x320(%r8),%xmm14
   0x00000000000013e1 <+5089>:	mulps  0x3a0(%rcx),%xmm14
   0x00000000000013e9 <+5097>:	addps  %xmm12,%xmm14
   0x0000000000001429 <+5161>:	movaps 0x360(%r8),%xmm14
   0x0000000000001431 <+5169>:	mulps  0x3a0(%rcx),%xmm14
   0x0000000000001439 <+5177>:	addps  %xmm15,%xmm14
   0x0000000000001479 <+5241>:	movaps 0x3a0(%r8),%xmm14
   0x0000000000001481 <+5249>:	mulps  0x3a0(%rcx),%xmm14
   0x0000000000001489 <+5257>:	addps  %xmm11,%xmm14
   0x00000000000014c9 <+5321>:	movaps 0x3e0(%r8),%xmm14
   0x00000000000014d1 <+5329>:	mulps  0x3a0(%rcx),%xmm14
   0x00000000000014d9 <+5337>:	addps  %xmm13,%xmm14

457
458	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
459	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
460	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000013ed <+5101>:	movaps 0x330(%r8),%xmm12
   0x00000000000013f5 <+5109>:	mulps  0x3b0(%rcx),%xmm12
   0x00000000000013fd <+5117>:	addps  %xmm14,%xmm12
   0x000000000000143d <+5181>:	movaps 0x370(%r8),%xmm15
   0x0000000000001445 <+5189>:	mulps  0x3b0(%rcx),%xmm15
   0x000000000000144d <+5197>:	addps  %xmm14,%xmm15
   0x000000000000148d <+5261>:	movaps 0x3b0(%r8),%xmm11
   0x0000000000001495 <+5269>:	mulps  0x3b0(%rcx),%xmm11
   0x000000000000149d <+5277>:	addps  %xmm14,%xmm11
   0x00000000000014dd <+5341>:	movaps 0x3f0(%r8),%xmm13
   0x00000000000014e5 <+5349>:	mulps  0x3b0(%rcx),%xmm13
   0x00000000000014ed <+5357>:	addps  %xmm14,%xmm13

461	      }
462	    }
463
464	    /* Store C(1,3) block. */
465	    for (i = 0; i < 4; i++)
466	    {
467	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000014f1 <+5361>:	mulps  %xmm0,%xmm12
   0x000000000000150d <+5389>:	mulps  %xmm0,%xmm15
   0x0000000000001521 <+5409>:	mulps  %xmm0,%xmm11
   0x0000000000001539 <+5433>:	mulps  %xmm0,%xmm13

468	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_13]), C_row[i]);
   0x00000000000014f5 <+5365>:	addps  0x80(%r9),%xmm12
   0x0000000000001511 <+5393>:	addps  0x90(%r9),%xmm15
   0x0000000000001525 <+5413>:	addps  0xa0(%r9),%xmm11
   0x000000000000153d <+5437>:	addps  0xb0(%r9),%xmm13

469	      _mm_store_ps(&C[i*4+C_OFFSET_13], C_row[i]);
   0x0000000000001501 <+5377>:	movaps %xmm12,0x80(%r9)
   0x0000000000001519 <+5401>:	movaps %xmm15,0x90(%r9)
   0x000000000000152d <+5421>:	movaps %xmm11,0xa0(%r9)
   0x0000000000001545 <+5445>:	movaps %xmm13,0xb0(%r9)

470	    }
471
472	    /* Reset C(1,4) matrix accumulators */
473	    C_row[0] = _mm_setzero_ps();
   0x0000000000001509 <+5385>:	xorps  %xmm12,%xmm12

474	    C_row[1] = _mm_setzero_ps();
   0x0000000000001535 <+5429>:	xorps  %xmm11,%xmm11

475	    C_row[2] = _mm_setzero_ps();
   0x00000000000014fd <+5373>:	xorps  %xmm14,%xmm14

476	    C_row[3] = _mm_setzero_ps();
   0x000000000000154d <+5453>:	xorps  %xmm13,%xmm13

477
478	    if (norm_product[0][3] &&
   0x0000000000001551 <+5457>:	comiss %xmm1,%xmm10
   0x0000000000001555 <+5461>:	jb     0x1a6e <stream_kernel+6766>

479	        norm_product[1][7] &&
   0x000000000000155b <+5467>:	movss  0x140(%rsp),%xmm10
   0x0000000000001565 <+5477>:	comiss %xmm1,%xmm10
   0x0000000000001569 <+5481>:	jb     0x1a6e <stream_kernel+6766>

480	        norm_product[2][11] &&
   0x000000000000156f <+5487>:	movss  0x138(%rsp),%xmm10
   0x0000000000001579 <+5497>:	comiss %xmm1,%xmm10
   0x000000000000157d <+5501>:	jb     0x1a6e <stream_kernel+6766>

481	        norm_product[3][15])
   0x0000000000001583 <+5507>:	movss  0x130(%rsp),%xmm10
   0x000000000000158d <+5517>:	comiss %xmm1,%xmm10
   0x0000000000001591 <+5521>:	jb     0x1a6e <stream_kernel+6766>

482	    {
483	      /* A(1,1)*B(1,4) = C(1,4). */
484	      for (i = 0; i < 4; i++)
485	      {
486	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
487	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
488	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001597 <+5527>:	movaps (%r8),%xmm12
   0x00000000000015aa <+5546>:	movaps 0x40(%r8),%xmm13
   0x00000000000015b4 <+5556>:	mulps  0xc0(%rcx),%xmm12
   0x00000000000015d4 <+5588>:	mulps  0xc0(%rcx),%xmm13
   0x000000000000160a <+5642>:	movaps 0x80(%r8),%xmm14
   0x0000000000001612 <+5650>:	mulps  0xc0(%rcx),%xmm14
   0x000000000000166e <+5742>:	movaps 0xc0(%r8),%xmm12
   0x0000000000001676 <+5750>:	mulps  0xc0(%rcx),%xmm12

489
490	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
491	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
492	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000159b <+5531>:	movaps 0x10(%r8),%xmm11
   0x00000000000015bc <+5564>:	mulps  0xd0(%rcx),%xmm11
   0x00000000000015e4 <+5604>:	addps  %xmm12,%xmm11
   0x00000000000015e8 <+5608>:	movaps 0x50(%r8),%xmm12
   0x00000000000015ed <+5613>:	mulps  0xd0(%rcx),%xmm12
   0x000000000000161a <+5658>:	addps  %xmm13,%xmm12
   0x0000000000001632 <+5682>:	movaps 0x90(%r8),%xmm12
   0x000000000000163a <+5690>:	mulps  0xd0(%rcx),%xmm12
   0x0000000000001656 <+5718>:	addps  %xmm14,%xmm12
   0x000000000000165a <+5722>:	movaps 0xd0(%r8),%xmm14
   0x0000000000001662 <+5730>:	mulps  0xd0(%rcx),%xmm14
   0x0000000000001692 <+5778>:	addps  %xmm12,%xmm14

493
494	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
495	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
496	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000015a0 <+5536>:	movaps 0x20(%r8),%xmm14
   0x00000000000015c4 <+5572>:	mulps  0xe0(%rcx),%xmm14
   0x00000000000015f5 <+5621>:	addps  %xmm11,%xmm14
   0x00000000000015f9 <+5625>:	movaps 0x60(%r8),%xmm11
   0x00000000000015fe <+5630>:	mulps  0xe0(%rcx),%xmm11
   0x000000000000161e <+5662>:	movaps 0xa0(%r8),%xmm13
   0x0000000000001626 <+5670>:	mulps  0xe0(%rcx),%xmm13
   0x000000000000162e <+5678>:	addps  %xmm12,%xmm11
   0x000000000000166a <+5738>:	addps  %xmm12,%xmm13
   0x0000000000001682 <+5762>:	movaps 0xe0(%r8),%xmm13
   0x000000000000168a <+5770>:	mulps  0xe0(%rcx),%xmm13
   0x00000000000016a6 <+5798>:	addps  %xmm14,%xmm13

497
498	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
499	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
500	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000015a5 <+5541>:	movaps 0x30(%r8),%xmm15
   0x00000000000015af <+5551>:	movaps 0x70(%r8),%xmm10
   0x00000000000015cc <+5580>:	mulps  0xf0(%rcx),%xmm15
   0x00000000000015dc <+5596>:	mulps  0xf0(%rcx),%xmm10
   0x0000000000001606 <+5638>:	addps  %xmm14,%xmm15
   0x0000000000001642 <+5698>:	addps  %xmm11,%xmm10
   0x0000000000001646 <+5702>:	movaps 0xb0(%r8),%xmm11
   0x000000000000164e <+5710>:	mulps  0xf0(%rcx),%xmm11
   0x000000000000167e <+5758>:	addps  %xmm13,%xmm11
   0x0000000000001696 <+5782>:	movaps 0xf0(%r8),%xmm12
   0x000000000000169e <+5790>:	mulps  0xf0(%rcx),%xmm12
   0x00000000000016ba <+5818>:	addps  %xmm13,%xmm12

501	      }
502
503	      /* A(1,2)*B(2,4) = C(1,4). */
504	      for (i = 0; i < 4; i++)
505	      {
506	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
507	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
508	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000016aa <+5802>:	movaps 0x100(%r8),%xmm14
   0x00000000000016b2 <+5810>:	mulps  0x1c0(%rcx),%xmm14
   0x00000000000016ce <+5838>:	addps  %xmm15,%xmm14
   0x00000000000016fa <+5882>:	movaps 0x140(%r8),%xmm15
   0x0000000000001702 <+5890>:	mulps  0x1c0(%rcx),%xmm15
   0x000000000000171e <+5918>:	addps  %xmm10,%xmm15
   0x000000000000174a <+5962>:	movaps 0x180(%r8),%xmm10
   0x0000000000001752 <+5970>:	mulps  0x1c0(%rcx),%xmm10
   0x000000000000176e <+5998>:	addps  %xmm11,%xmm10
   0x000000000000179a <+6042>:	movaps 0x1c0(%r8),%xmm11
   0x00000000000017a2 <+6050>:	mulps  0x1c0(%rcx),%xmm11
   0x00000000000017be <+6078>:	addps  %xmm12,%xmm11

509
510	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
511	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
512	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000016d2 <+5842>:	movaps 0x110(%r8),%xmm15
   0x00000000000016da <+5850>:	mulps  0x1d0(%rcx),%xmm15
   0x00000000000016e2 <+5858>:	addps  %xmm14,%xmm15
   0x0000000000001722 <+5922>:	movaps 0x150(%r8),%xmm10
   0x000000000000172a <+5930>:	mulps  0x1d0(%rcx),%xmm10
   0x0000000000001732 <+5938>:	addps  %xmm15,%xmm10
   0x0000000000001772 <+6002>:	movaps 0x190(%r8),%xmm11
   0x000000000000177a <+6010>:	mulps  0x1d0(%rcx),%xmm11
   0x0000000000001782 <+6018>:	addps  %xmm10,%xmm11
   0x00000000000017c2 <+6082>:	movaps 0x1d0(%r8),%xmm12
   0x00000000000017ca <+6090>:	mulps  0x1d0(%rcx),%xmm12
   0x00000000000017d2 <+6098>:	addps  %xmm11,%xmm12

513
514	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
515	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
516	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000016be <+5822>:	movaps 0x120(%r8),%xmm13
   0x00000000000016c6 <+5830>:	mulps  0x1e0(%rcx),%xmm13
   0x00000000000016f6 <+5878>:	addps  %xmm15,%xmm13
   0x0000000000001736 <+5942>:	movaps 0x160(%r8),%xmm15
   0x000000000000173e <+5950>:	mulps  0x1e0(%rcx),%xmm15
   0x0000000000001746 <+5958>:	addps  %xmm10,%xmm15
   0x000000000000175e <+5982>:	movaps 0x1a0(%r8),%xmm15
   0x0000000000001766 <+5990>:	mulps  0x1e0(%rcx),%xmm15
   0x0000000000001796 <+6038>:	addps  %xmm11,%xmm15
   0x00000000000017ae <+6062>:	movaps 0x1e0(%r8),%xmm15
   0x00000000000017b6 <+6070>:	mulps  0x1e0(%rcx),%xmm15
   0x00000000000017e6 <+6118>:	addps  %xmm12,%xmm15

517
518	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
519	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
520	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000016e6 <+5862>:	movaps 0x130(%r8),%xmm14
   0x00000000000016ee <+5870>:	mulps  0x1f0(%rcx),%xmm14
   0x000000000000170a <+5898>:	addps  %xmm13,%xmm14
   0x000000000000170e <+5902>:	movaps 0x170(%r8),%xmm13
   0x0000000000001716 <+5910>:	mulps  0x1f0(%rcx),%xmm13
   0x000000000000175a <+5978>:	addps  %xmm15,%xmm13
   0x0000000000001786 <+6022>:	movaps 0x1b0(%r8),%xmm10
   0x000000000000178e <+6030>:	mulps  0x1f0(%rcx),%xmm10
   0x00000000000017aa <+6058>:	addps  %xmm15,%xmm10
   0x00000000000017d6 <+6102>:	movaps 0x1f0(%r8),%xmm11
   0x00000000000017de <+6110>:	mulps  0x1f0(%rcx),%xmm11
   0x00000000000017fa <+6138>:	addps  %xmm15,%xmm11

521	      }
522
523	      /* A(1,3)*B(3,4) = C(1,4). */
524	      for (i = 0; i < 4; i++)
525	      {
526	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
527	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
528	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000017ea <+6122>:	movaps 0x200(%r8),%xmm12
   0x00000000000017f2 <+6130>:	mulps  0x2c0(%rcx),%xmm12
   0x000000000000180e <+6158>:	addps  %xmm14,%xmm12
   0x000000000000183a <+6202>:	movaps 0x240(%r8),%xmm14
   0x0000000000001842 <+6210>:	mulps  0x2c0(%rcx),%xmm14
   0x000000000000185e <+6238>:	addps  %xmm13,%xmm14
   0x0000000000001876 <+6262>:	movaps 0x280(%r8),%xmm14
   0x000000000000187e <+6270>:	mulps  0x2c0(%rcx),%xmm14
   0x000000000000189a <+6298>:	addps  %xmm10,%xmm14
   0x00000000000018da <+6362>:	movaps 0x2c0(%r8),%xmm10
   0x00000000000018e2 <+6370>:	mulps  0x2c0(%rcx),%xmm10
   0x00000000000018fe <+6398>:	addps  %xmm11,%xmm10

529
530	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
531	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
532	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001812 <+6162>:	movaps 0x210(%r8),%xmm14
   0x000000000000181a <+6170>:	mulps  0x2d0(%rcx),%xmm14
   0x0000000000001822 <+6178>:	addps  %xmm12,%xmm14
   0x0000000000001862 <+6242>:	movaps 0x250(%r8),%xmm13
   0x000000000000186a <+6250>:	mulps  0x2d0(%rcx),%xmm13
   0x0000000000001872 <+6258>:	addps  %xmm14,%xmm13
   0x000000000000189e <+6302>:	movaps 0x290(%r8),%xmm10
   0x00000000000018a6 <+6310>:	mulps  0x2d0(%rcx),%xmm10
   0x00000000000018c2 <+6338>:	addps  %xmm14,%xmm10
   0x0000000000001902 <+6402>:	movaps 0x2d0(%r8),%xmm11
   0x000000000000190a <+6410>:	mulps  0x2d0(%rcx),%xmm11
   0x0000000000001912 <+6418>:	addps  %xmm10,%xmm11

533
534	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
535	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
536	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000017fe <+6142>:	movaps 0x220(%r8),%xmm15
   0x0000000000001806 <+6150>:	mulps  0x2e0(%rcx),%xmm15
   0x0000000000001836 <+6198>:	addps  %xmm14,%xmm15
   0x000000000000184e <+6222>:	movaps 0x260(%r8),%xmm15
   0x0000000000001856 <+6230>:	mulps  0x2e0(%rcx),%xmm15
   0x0000000000001886 <+6278>:	addps  %xmm13,%xmm15
   0x00000000000018b2 <+6322>:	movaps 0x2a0(%r8),%xmm15
   0x00000000000018ba <+6330>:	mulps  0x2e0(%rcx),%xmm15
   0x00000000000018d6 <+6358>:	addps  %xmm10,%xmm15
   0x0000000000001916 <+6422>:	movaps 0x2e0(%r8),%xmm10
   0x000000000000191e <+6430>:	mulps  0x2e0(%rcx),%xmm10
   0x0000000000001926 <+6438>:	addps  %xmm11,%xmm10

537
538	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
539	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
540	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001826 <+6182>:	movaps 0x230(%r8),%xmm12
   0x000000000000182e <+6190>:	mulps  0x2f0(%rcx),%xmm12
   0x000000000000184a <+6218>:	addps  %xmm15,%xmm12
   0x000000000000188a <+6282>:	movaps 0x270(%r8),%xmm13
   0x0000000000001892 <+6290>:	mulps  0x2f0(%rcx),%xmm13
   0x00000000000018ae <+6318>:	addps  %xmm15,%xmm13
   0x00000000000018c6 <+6342>:	movaps 0x2b0(%r8),%xmm14
   0x00000000000018ce <+6350>:	mulps  0x2f0(%rcx),%xmm14
   0x00000000000018ea <+6378>:	addps  %xmm15,%xmm14
   0x00000000000018ee <+6382>:	movaps 0x2f0(%r8),%xmm15
   0x00000000000018f6 <+6390>:	mulps  0x2f0(%rcx),%xmm15
   0x000000000000193a <+6458>:	addps  %xmm10,%xmm15

541	      }
542
543	      /* A(1,4)*B(4,4) = C(1,4). */
544	      for (i = 0; i < 4; i++)
545	      {
546	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
547	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
548	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000192a <+6442>:	movaps 0x300(%r8),%xmm11
   0x0000000000001932 <+6450>:	mulps  0x3c0(%rcx),%xmm11
   0x000000000000194e <+6478>:	addps  %xmm12,%xmm11
   0x0000000000001966 <+6502>:	movaps 0x340(%r8),%xmm11
   0x000000000000196e <+6510>:	mulps  0x3c0(%rcx),%xmm11
   0x000000000000198a <+6538>:	addps  %xmm13,%xmm11
   0x00000000000019ca <+6602>:	movaps 0x380(%r8),%xmm13
   0x00000000000019d2 <+6610>:	mulps  0x3c0(%rcx),%xmm13
   0x00000000000019ee <+6638>:	addps  %xmm14,%xmm13
   0x0000000000001a06 <+6662>:	movaps 0x3c0(%r8),%xmm13
   0x0000000000001a0e <+6670>:	mulps  0x3c0(%rcx),%xmm13
   0x0000000000001a2a <+6698>:	addps  %xmm15,%xmm13

549
550	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
551	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
552	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001952 <+6482>:	movaps 0x310(%r8),%xmm12
   0x000000000000195a <+6490>:	mulps  0x3d0(%rcx),%xmm12
   0x0000000000001962 <+6498>:	addps  %xmm11,%xmm12
   0x000000000000198e <+6542>:	movaps 0x350(%r8),%xmm13
   0x0000000000001996 <+6550>:	mulps  0x3d0(%rcx),%xmm13
   0x00000000000019b2 <+6578>:	addps  %xmm11,%xmm13
   0x00000000000019f2 <+6642>:	movaps 0x390(%r8),%xmm14
   0x00000000000019fa <+6650>:	mulps  0x3d0(%rcx),%xmm14
   0x0000000000001a02 <+6658>:	addps  %xmm13,%xmm14
   0x0000000000001a2e <+6702>:	movaps 0x3d0(%r8),%xmm15
   0x0000000000001a36 <+6710>:	mulps  0x3d0(%rcx),%xmm15
   0x0000000000001a52 <+6738>:	addps  %xmm13,%xmm15

553
554	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
555	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
556	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000193e <+6462>:	movaps 0x320(%r8),%xmm10
   0x0000000000001946 <+6470>:	mulps  0x3e0(%rcx),%xmm10
   0x0000000000001976 <+6518>:	addps  %xmm12,%xmm10
   0x00000000000019a2 <+6562>:	movaps 0x360(%r8),%xmm10
   0x00000000000019aa <+6570>:	mulps  0x3e0(%rcx),%xmm10
   0x00000000000019c6 <+6598>:	addps  %xmm13,%xmm10
   0x00000000000019de <+6622>:	movaps 0x3a0(%r8),%xmm10
   0x00000000000019e6 <+6630>:	mulps  0x3e0(%rcx),%xmm10
   0x0000000000001a16 <+6678>:	addps  %xmm14,%xmm10
   0x0000000000001a42 <+6722>:	movaps 0x3e0(%r8),%xmm10
   0x0000000000001a4a <+6730>:	mulps  0x3e0(%rcx),%xmm10
   0x0000000000001a66 <+6758>:	addps  %xmm15,%xmm10

557
558	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
559	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
560	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000197a <+6522>:	movaps 0x330(%r8),%xmm12
   0x0000000000001982 <+6530>:	mulps  0x3f0(%rcx),%xmm12
   0x000000000000199e <+6558>:	addps  %xmm10,%xmm12
   0x00000000000019b6 <+6582>:	movaps 0x370(%r8),%xmm11
   0x00000000000019be <+6590>:	mulps  0x3f0(%rcx),%xmm11
   0x00000000000019da <+6618>:	addps  %xmm10,%xmm11
   0x0000000000001a1a <+6682>:	movaps 0x3b0(%r8),%xmm14
   0x0000000000001a22 <+6690>:	mulps  0x3f0(%rcx),%xmm14
   0x0000000000001a3e <+6718>:	addps  %xmm10,%xmm14
   0x0000000000001a56 <+6742>:	movaps 0x3f0(%r8),%xmm13
   0x0000000000001a5e <+6750>:	mulps  0x3f0(%rcx),%xmm13
   0x0000000000001a6a <+6762>:	addps  %xmm10,%xmm13

561	      }
562	    }
563
564	    /* Store C(1,4) block. */
565	    for (i = 0; i < 4; i++)
566	    {
567	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000001a6e <+6766>:	mulps  %xmm0,%xmm12
   0x0000000000001a8a <+6794>:	mulps  %xmm0,%xmm11
   0x0000000000001aa2 <+6818>:	mulps  %xmm0,%xmm14
   0x0000000000001ac0 <+6848>:	mulps  %xmm0,%xmm13

568	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_14]), C_row[i]);
   0x0000000000001a72 <+6770>:	addps  0xc0(%r9),%xmm12
   0x0000000000001a8e <+6798>:	addps  0xd0(%r9),%xmm11
   0x0000000000001aa6 <+6822>:	addps  0xe0(%r9),%xmm14
   0x0000000000001ac4 <+6852>:	addps  0xf0(%r9),%xmm13

569	      _mm_store_ps(&C[i*4+C_OFFSET_14], C_row[i]);
   0x0000000000001a7e <+6782>:	movaps %xmm12,0xc0(%r9)
   0x0000000000001a96 <+6806>:	movaps %xmm11,0xd0(%r9)
   0x0000000000001aae <+6830>:	movaps %xmm14,0xe0(%r9)
   0x0000000000001acc <+6860>:	movaps %xmm13,0xf0(%r9)

570	    }
571
572	    /* Reset C(2,1) matrix accumulators */
573	    C_row[0] = _mm_setzero_ps();
   0x0000000000001a7a <+6778>:	xorps  %xmm10,%xmm10

574	    C_row[1] = _mm_setzero_ps();
   0x0000000000001ad4 <+6868>:	xorps  %xmm13,%xmm13

575	    C_row[2] = _mm_setzero_ps();
   0x0000000000001a86 <+6790>:	xorps  %xmm12,%xmm12

576	    C_row[3] = _mm_setzero_ps();
   0x0000000000001a9e <+6814>:	xorps  %xmm11,%xmm11

577
578	    if (norm_product[4][0] &&
   0x0000000000001ab6 <+6838>:	movss  0x170(%rsp),%xmm14
   0x0000000000001ad8 <+6872>:	comiss %xmm1,%xmm14
   0x0000000000001adc <+6876>:	jb     0x1fda <stream_kernel+8154>

579	        norm_product[5][4] &&
   0x0000000000001ae2 <+6882>:	movss  0x120(%rsp),%xmm14
   0x0000000000001aec <+6892>:	comiss %xmm1,%xmm14
   0x0000000000001af0 <+6896>:	jb     0x1fda <stream_kernel+8154>

580	        norm_product[6][8] &&
   0x0000000000001af6 <+6902>:	movss  0x118(%rsp),%xmm14
   0x0000000000001b00 <+6912>:	comiss %xmm1,%xmm14
   0x0000000000001b04 <+6916>:	jb     0x1fda <stream_kernel+8154>

581	        norm_product[7][12])
   0x0000000000001b0a <+6922>:	movss  0x110(%rsp),%xmm14
   0x0000000000001b14 <+6932>:	comiss %xmm1,%xmm14
   0x0000000000001b18 <+6936>:	jb     0x1fda <stream_kernel+8154>

582	    {
583	      /* A(2,1)*B(1,1) = C(2,1). */
584	      for (i = 0; i < 4; i++)
585	      {
586	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
587	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
588	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b1e <+6942>:	movaps 0x400(%r8),%xmm10
   0x0000000000001b3e <+6974>:	movaps 0x440(%r8),%xmm15
   0x0000000000001b4e <+6990>:	mulps  (%rcx),%xmm10
   0x0000000000001b61 <+7009>:	mulps  (%rcx),%xmm15
   0x0000000000001b7f <+7039>:	movaps 0x480(%r8),%xmm13
   0x0000000000001b87 <+7047>:	mulps  (%rcx),%xmm13
   0x0000000000001bb1 <+7089>:	movaps 0x4c0(%r8),%xmm14
   0x0000000000001bb9 <+7097>:	mulps  (%rcx),%xmm14

589
590	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
591	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
592	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b26 <+6950>:	movaps 0x410(%r8),%xmm13
   0x0000000000001b46 <+6982>:	movaps 0x450(%r8),%xmm14
   0x0000000000001b52 <+6994>:	mulps  0x10(%rcx),%xmm13
   0x0000000000001b65 <+7013>:	mulps  0x10(%rcx),%xmm14
   0x0000000000001b6a <+7018>:	addps  %xmm10,%xmm13
   0x0000000000001b9c <+7068>:	addps  %xmm15,%xmm14
   0x0000000000001ba0 <+7072>:	movaps 0x490(%r8),%xmm15
   0x0000000000001ba8 <+7080>:	mulps  0x10(%rcx),%xmm15
   0x0000000000001bce <+7118>:	addps  %xmm13,%xmm15
   0x0000000000001bf4 <+7156>:	movaps 0x4d0(%r8),%xmm13
   0x0000000000001bfc <+7164>:	mulps  0x10(%rcx),%xmm13
   0x0000000000001c01 <+7169>:	addps  %xmm14,%xmm13

593
594	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
595	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
596	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b2e <+6958>:	movaps 0x420(%r8),%xmm12
   0x0000000000001b57 <+6999>:	mulps  0x20(%rcx),%xmm12
   0x0000000000001b6e <+7022>:	movaps 0x460(%r8),%xmm10
   0x0000000000001b76 <+7030>:	mulps  0x20(%rcx),%xmm10
   0x0000000000001b7b <+7035>:	addps  %xmm13,%xmm12
   0x0000000000001bad <+7085>:	addps  %xmm14,%xmm10
   0x0000000000001bd2 <+7122>:	movaps 0x4a0(%r8),%xmm13
   0x0000000000001bda <+7130>:	mulps  0x20(%rcx),%xmm13
   0x0000000000001bdf <+7135>:	addps  %xmm15,%xmm13
   0x0000000000001be3 <+7139>:	movaps 0x4e0(%r8),%xmm15
   0x0000000000001beb <+7147>:	mulps  0x20(%rcx),%xmm15
   0x0000000000001c12 <+7186>:	addps  %xmm13,%xmm15

597
598	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
599	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
600	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b36 <+6966>:	movaps 0x430(%r8),%xmm11
   0x0000000000001b5c <+7004>:	mulps  0x30(%rcx),%xmm11
   0x0000000000001b8b <+7051>:	addps  %xmm12,%xmm11
   0x0000000000001b8f <+7055>:	movaps 0x470(%r8),%xmm12
   0x0000000000001b97 <+7063>:	mulps  0x30(%rcx),%xmm12
   0x0000000000001bbd <+7101>:	addps  %xmm10,%xmm12
   0x0000000000001bc1 <+7105>:	movaps 0x4b0(%r8),%xmm10
   0x0000000000001bc9 <+7113>:	mulps  0x30(%rcx),%xmm10
   0x0000000000001bf0 <+7152>:	addps  %xmm13,%xmm10
   0x0000000000001c05 <+7173>:	movaps 0x4f0(%r8),%xmm14
   0x0000000000001c0d <+7181>:	mulps  0x30(%rcx),%xmm14
   0x0000000000001c26 <+7206>:	addps  %xmm15,%xmm14

601	      }
602
603	      /* A(2,2)*B(2,1) = C(2,1). */
604	      for (i = 0; i < 4; i++)
605	      {
606	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
607	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
608	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c16 <+7190>:	movaps 0x500(%r8),%xmm13
   0x0000000000001c1e <+7198>:	mulps  0x100(%rcx),%xmm13
   0x0000000000001c3a <+7226>:	addps  %xmm11,%xmm13
   0x0000000000001c66 <+7270>:	movaps 0x540(%r8),%xmm11
   0x0000000000001c6e <+7278>:	mulps  0x100(%rcx),%xmm11
   0x0000000000001c8a <+7306>:	addps  %xmm12,%xmm11
   0x0000000000001cb6 <+7350>:	movaps 0x580(%r8),%xmm12
   0x0000000000001cbe <+7358>:	mulps  0x100(%rcx),%xmm12
   0x0000000000001cda <+7386>:	addps  %xmm10,%xmm12
   0x0000000000001cf2 <+7410>:	movaps 0x5c0(%r8),%xmm12
   0x0000000000001cfa <+7418>:	mulps  0x100(%rcx),%xmm12
   0x0000000000001d16 <+7446>:	addps  %xmm14,%xmm12

609
610	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
611	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
612	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c3e <+7230>:	movaps 0x510(%r8),%xmm11
   0x0000000000001c46 <+7238>:	mulps  0x110(%rcx),%xmm11
   0x0000000000001c4e <+7246>:	addps  %xmm13,%xmm11
   0x0000000000001c8e <+7310>:	movaps 0x550(%r8),%xmm12
   0x0000000000001c96 <+7318>:	mulps  0x110(%rcx),%xmm12
   0x0000000000001c9e <+7326>:	addps  %xmm11,%xmm12
   0x0000000000001cde <+7390>:	movaps 0x590(%r8),%xmm10
   0x0000000000001ce6 <+7398>:	mulps  0x110(%rcx),%xmm10
   0x0000000000001cee <+7406>:	addps  %xmm12,%xmm10
   0x0000000000001d1a <+7450>:	movaps 0x5d0(%r8),%xmm14
   0x0000000000001d22 <+7458>:	mulps  0x110(%rcx),%xmm14
   0x0000000000001d3e <+7486>:	addps  %xmm12,%xmm14

613
614	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
615	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
616	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c2a <+7210>:	movaps 0x520(%r8),%xmm15
   0x0000000000001c32 <+7218>:	mulps  0x120(%rcx),%xmm15
   0x0000000000001c62 <+7266>:	addps  %xmm11,%xmm15
   0x0000000000001ca2 <+7330>:	movaps 0x560(%r8),%xmm11
   0x0000000000001caa <+7338>:	mulps  0x120(%rcx),%xmm11
   0x0000000000001cb2 <+7346>:	addps  %xmm12,%xmm11
   0x0000000000001cca <+7370>:	movaps 0x5a0(%r8),%xmm11
   0x0000000000001cd2 <+7378>:	mulps  0x120(%rcx),%xmm11
   0x0000000000001d02 <+7426>:	addps  %xmm10,%xmm11
   0x0000000000001d42 <+7490>:	movaps 0x5e0(%r8),%xmm12
   0x0000000000001d4a <+7498>:	mulps  0x120(%rcx),%xmm12
   0x0000000000001d52 <+7506>:	addps  %xmm14,%xmm12

617
618	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
619	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
620	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c52 <+7250>:	movaps 0x530(%r8),%xmm13
   0x0000000000001c5a <+7258>:	mulps  0x130(%rcx),%xmm13
   0x0000000000001c76 <+7286>:	addps  %xmm15,%xmm13
   0x0000000000001c7a <+7290>:	movaps 0x570(%r8),%xmm15
   0x0000000000001c82 <+7298>:	mulps  0x130(%rcx),%xmm15
   0x0000000000001cc6 <+7366>:	addps  %xmm11,%xmm15
   0x0000000000001d06 <+7430>:	movaps 0x5b0(%r8),%xmm10
   0x0000000000001d0e <+7438>:	mulps  0x130(%rcx),%xmm10
   0x0000000000001d2a <+7466>:	addps  %xmm11,%xmm10
   0x0000000000001d2e <+7470>:	movaps 0x5f0(%r8),%xmm11
   0x0000000000001d36 <+7478>:	mulps  0x130(%rcx),%xmm11
   0x0000000000001d66 <+7526>:	addps  %xmm12,%xmm11

621	      }
622
623	      /* A(2,3)*B(3,1) = C(2,1). */
624	      for (i = 0; i < 4; i++)
625	      {
626	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
627	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
628	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001d56 <+7510>:	movaps 0x600(%r8),%xmm14
   0x0000000000001d5e <+7518>:	mulps  0x200(%rcx),%xmm14
   0x0000000000001d7a <+7546>:	addps  %xmm13,%xmm14
   0x0000000000001da6 <+7590>:	movaps 0x640(%r8),%xmm13
   0x0000000000001dae <+7598>:	mulps  0x200(%rcx),%xmm13
   0x0000000000001dca <+7626>:	addps  %xmm15,%xmm13
   0x0000000000001df6 <+7670>:	movaps 0x680(%r8),%xmm15
   0x0000000000001dfe <+7678>:	mulps  0x200(%rcx),%xmm15
   0x0000000000001e1a <+7706>:	addps  %xmm10,%xmm15
   0x0000000000001e46 <+7750>:	movaps 0x6c0(%r8),%xmm10
   0x0000000000001e4e <+7758>:	mulps  0x200(%rcx),%xmm10
   0x0000000000001e6a <+7786>:	addps  %xmm11,%xmm10

629
630	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
631	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
632	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001d7e <+7550>:	movaps 0x610(%r8),%xmm13
   0x0000000000001d86 <+7558>:	mulps  0x210(%rcx),%xmm13
   0x0000000000001d8e <+7566>:	addps  %xmm14,%xmm13
   0x0000000000001dce <+7630>:	movaps 0x650(%r8),%xmm15
   0x0000000000001dd6 <+7638>:	mulps  0x210(%rcx),%xmm15
   0x0000000000001dde <+7646>:	addps  %xmm13,%xmm15
   0x0000000000001e1e <+7710>:	movaps 0x690(%r8),%xmm10
   0x0000000000001e26 <+7718>:	mulps  0x210(%rcx),%xmm10
   0x0000000000001e2e <+7726>:	addps  %xmm15,%xmm10
   0x0000000000001e6e <+7790>:	movaps 0x6d0(%r8),%xmm11
   0x0000000000001e76 <+7798>:	mulps  0x210(%rcx),%xmm11
   0x0000000000001e7e <+7806>:	addps  %xmm10,%xmm11

633
634	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
635	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
636	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001d92 <+7570>:	movaps 0x620(%r8),%xmm14
   0x0000000000001d9a <+7578>:	mulps  0x220(%rcx),%xmm14
   0x0000000000001da2 <+7586>:	addps  %xmm13,%xmm14
   0x0000000000001dba <+7610>:	movaps 0x660(%r8),%xmm14
   0x0000000000001dc2 <+7618>:	mulps  0x220(%rcx),%xmm14
   0x0000000000001df2 <+7666>:	addps  %xmm15,%xmm14
   0x0000000000001e32 <+7730>:	movaps 0x6a0(%r8),%xmm15
   0x0000000000001e3a <+7738>:	mulps  0x220(%rcx),%xmm15
   0x0000000000001e42 <+7746>:	addps  %xmm10,%xmm15
   0x0000000000001e82 <+7810>:	movaps 0x6e0(%r8),%xmm10
   0x0000000000001e8a <+7818>:	mulps  0x220(%rcx),%xmm10
   0x0000000000001e92 <+7826>:	addps  %xmm11,%xmm10

637
638	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
639	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
640	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001d6a <+7530>:	movaps 0x630(%r8),%xmm12
   0x0000000000001d72 <+7538>:	mulps  0x230(%rcx),%xmm12
   0x0000000000001db6 <+7606>:	addps  %xmm14,%xmm12
   0x0000000000001de2 <+7650>:	movaps 0x670(%r8),%xmm13
   0x0000000000001dea <+7658>:	mulps  0x230(%rcx),%xmm13
   0x0000000000001e06 <+7686>:	addps  %xmm14,%xmm13
   0x0000000000001e0a <+7690>:	movaps 0x6b0(%r8),%xmm14
   0x0000000000001e12 <+7698>:	mulps  0x230(%rcx),%xmm14
   0x0000000000001e56 <+7766>:	addps  %xmm15,%xmm14
   0x0000000000001e5a <+7770>:	movaps 0x6f0(%r8),%xmm15
   0x0000000000001e62 <+7778>:	mulps  0x230(%rcx),%xmm15
   0x0000000000001ea6 <+7846>:	addps  %xmm10,%xmm15

641	      }
642
643	      /* A(2,4)*B(4,1) = C(2,1). */
644	      for (i = 0; i < 4; i++)
645	      {
646	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
647	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
648	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001e96 <+7830>:	movaps 0x700(%r8),%xmm11
   0x0000000000001e9e <+7838>:	mulps  0x300(%rcx),%xmm11
   0x0000000000001eba <+7866>:	addps  %xmm12,%xmm11
   0x0000000000001ee6 <+7910>:	movaps 0x740(%r8),%xmm12
   0x0000000000001eee <+7918>:	mulps  0x300(%rcx),%xmm12
   0x0000000000001f0a <+7946>:	addps  %xmm13,%xmm12
   0x0000000000001f22 <+7970>:	movaps 0x780(%r8),%xmm12
   0x0000000000001f2a <+7978>:	mulps  0x300(%rcx),%xmm12
   0x0000000000001f46 <+8006>:	addps  %xmm14,%xmm12
   0x0000000000001f86 <+8070>:	movaps 0x7c0(%r8),%xmm14
   0x0000000000001f8e <+8078>:	mulps  0x300(%rcx),%xmm14
   0x0000000000001faa <+8106>:	addps  %xmm15,%xmm14

649
650	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
651	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
652	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001ebe <+7870>:	movaps 0x710(%r8),%xmm12
   0x0000000000001ec6 <+7878>:	mulps  0x310(%rcx),%xmm12
   0x0000000000001ece <+7886>:	addps  %xmm11,%xmm12
   0x0000000000001f0e <+7950>:	movaps 0x750(%r8),%xmm13
   0x0000000000001f16 <+7958>:	mulps  0x310(%rcx),%xmm13
   0x0000000000001f1e <+7966>:	addps  %xmm12,%xmm13
   0x0000000000001f4a <+8010>:	movaps 0x790(%r8),%xmm14
   0x0000000000001f52 <+8018>:	mulps  0x310(%rcx),%xmm14
   0x0000000000001f6e <+8046>:	addps  %xmm12,%xmm14
   0x0000000000001fae <+8110>:	movaps 0x7d0(%r8),%xmm15
   0x0000000000001fb6 <+8118>:	mulps  0x310(%rcx),%xmm15
   0x0000000000001fbe <+8126>:	addps  %xmm14,%xmm15

653
654	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
655	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
656	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001ed2 <+7890>:	movaps 0x720(%r8),%xmm11
   0x0000000000001eda <+7898>:	mulps  0x320(%rcx),%xmm11
   0x0000000000001ee2 <+7906>:	addps  %xmm12,%xmm11
   0x0000000000001efa <+7930>:	movaps 0x760(%r8),%xmm11
   0x0000000000001f02 <+7938>:	mulps  0x320(%rcx),%xmm11
   0x0000000000001f32 <+7986>:	addps  %xmm13,%xmm11
   0x0000000000001f5e <+8030>:	movaps 0x7a0(%r8),%xmm11
   0x0000000000001f66 <+8038>:	mulps  0x320(%rcx),%xmm11
   0x0000000000001f82 <+8066>:	addps  %xmm14,%xmm11
   0x0000000000001fc2 <+8130>:	movaps 0x7e0(%r8),%xmm14
   0x0000000000001fca <+8138>:	mulps  0x320(%rcx),%xmm14
   0x0000000000001fd2 <+8146>:	addps  %xmm15,%xmm14

657
658	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
659	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
660	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001eaa <+7850>:	movaps 0x730(%r8),%xmm10
   0x0000000000001eb2 <+7858>:	mulps  0x330(%rcx),%xmm10
   0x0000000000001ef6 <+7926>:	addps  %xmm11,%xmm10
   0x0000000000001f36 <+7990>:	movaps 0x770(%r8),%xmm13
   0x0000000000001f3e <+7998>:	mulps  0x330(%rcx),%xmm13
   0x0000000000001f5a <+8026>:	addps  %xmm11,%xmm13
   0x0000000000001f72 <+8050>:	movaps 0x7b0(%r8),%xmm12
   0x0000000000001f7a <+8058>:	mulps  0x330(%rcx),%xmm12
   0x0000000000001f96 <+8086>:	addps  %xmm11,%xmm12
   0x0000000000001f9a <+8090>:	movaps 0x7f0(%r8),%xmm11
   0x0000000000001fa2 <+8098>:	mulps  0x330(%rcx),%xmm11
   0x0000000000001fd6 <+8150>:	addps  %xmm14,%xmm11

661	      }
662	    }
663
664	    /* Store C(2,1) block. */
665	    for (i = 0; i < 4; i++)
666	    {
667	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000001fda <+8154>:	mulps  %xmm0,%xmm10
   0x0000000000001ffc <+8188>:	mulps  %xmm0,%xmm13
   0x0000000000002014 <+8212>:	mulps  %xmm0,%xmm12
   0x000000000000202c <+8236>:	mulps  %xmm0,%xmm11

668	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_21]), C_row[i]);
   0x0000000000001fde <+8158>:	addps  0x100(%r9),%xmm10
   0x0000000000002000 <+8192>:	addps  0x110(%r9),%xmm13
   0x0000000000002018 <+8216>:	addps  0x120(%r9),%xmm12
   0x0000000000002030 <+8240>:	addps  0x130(%r9),%xmm11

669	      _mm_store_ps(&C[i*4+C_OFFSET_21], C_row[i]);
   0x0000000000001ff0 <+8176>:	movaps %xmm10,0x100(%r9)
   0x0000000000002008 <+8200>:	movaps %xmm13,0x110(%r9)
   0x0000000000002020 <+8224>:	movaps %xmm12,0x120(%r9)
   0x0000000000002038 <+8248>:	movaps %xmm11,0x130(%r9)

670	    }
671
672	    /* Reset C(2,2) matrix accumulators */
673	    C_row[0] = _mm_setzero_ps();
   0x0000000000001ff8 <+8184>:	xorps  %xmm10,%xmm10

674	    C_row[1] = _mm_setzero_ps();
   0x0000000000002010 <+8208>:	xorps  %xmm13,%xmm13

675	    C_row[2] = _mm_setzero_ps();
   0x0000000000002028 <+8232>:	xorps  %xmm12,%xmm12

676	    C_row[3] = _mm_setzero_ps();
   0x0000000000002040 <+8256>:	xorps  %xmm11,%xmm11

677
678	    if (norm_product[4][1] &&
   0x0000000000001fe6 <+8166>:	movss  0x160(%rsp),%xmm14
   0x0000000000002044 <+8260>:	comiss %xmm1,%xmm14
   0x0000000000002048 <+8264>:	jb     0x254a <stream_kernel+9546>

679	        norm_product[5][5] &&
   0x000000000000204e <+8270>:	movss  0x100(%rsp),%xmm14
   0x0000000000002058 <+8280>:	comiss %xmm1,%xmm14
   0x000000000000205c <+8284>:	jb     0x254a <stream_kernel+9546>

680	        norm_product[6][9] &&
   0x0000000000002062 <+8290>:	movss  0xf8(%rsp),%xmm14
   0x000000000000206c <+8300>:	comiss %xmm1,%xmm14
   0x0000000000002070 <+8304>:	jb     0x254a <stream_kernel+9546>

681	        norm_product[7][13])
   0x0000000000002076 <+8310>:	movss  0xf0(%rsp),%xmm14
   0x0000000000002080 <+8320>:	comiss %xmm1,%xmm14
   0x0000000000002084 <+8324>:	jb     0x254a <stream_kernel+9546>

682	    {
683	      /* A(2,1)*B(1,2) = C(2,2). */
684	      for (i = 0; i < 4; i++)
685	      {
686	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
687	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
688	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000208a <+8330>:	movaps 0x400(%r8),%xmm10
   0x00000000000020aa <+8362>:	movaps 0x440(%r8),%xmm15
   0x00000000000020ba <+8378>:	mulps  0x40(%rcx),%xmm10
   0x00000000000020ce <+8398>:	mulps  0x40(%rcx),%xmm15
   0x00000000000020ed <+8429>:	movaps 0x480(%r8),%xmm13
   0x00000000000020f5 <+8437>:	mulps  0x40(%rcx),%xmm13
   0x0000000000002120 <+8480>:	movaps 0x4c0(%r8),%xmm14
   0x0000000000002128 <+8488>:	mulps  0x40(%rcx),%xmm14

689
690	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
691	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
692	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002092 <+8338>:	movaps 0x410(%r8),%xmm13
   0x00000000000020b2 <+8370>:	movaps 0x450(%r8),%xmm14
   0x00000000000020bf <+8383>:	mulps  0x50(%rcx),%xmm13
   0x00000000000020d3 <+8403>:	mulps  0x50(%rcx),%xmm14
   0x00000000000020d8 <+8408>:	addps  %xmm10,%xmm13
   0x000000000000210b <+8459>:	addps  %xmm15,%xmm14
   0x000000000000210f <+8463>:	movaps 0x490(%r8),%xmm15
   0x0000000000002117 <+8471>:	mulps  0x50(%rcx),%xmm15
   0x000000000000213e <+8510>:	addps  %xmm13,%xmm15
   0x0000000000002164 <+8548>:	movaps 0x4d0(%r8),%xmm13
   0x000000000000216c <+8556>:	mulps  0x50(%rcx),%xmm13
   0x0000000000002171 <+8561>:	addps  %xmm14,%xmm13

693
694	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
695	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
696	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000209a <+8346>:	movaps 0x420(%r8),%xmm12
   0x00000000000020c4 <+8388>:	mulps  0x60(%rcx),%xmm12
   0x00000000000020dc <+8412>:	movaps 0x460(%r8),%xmm10
   0x00000000000020e4 <+8420>:	mulps  0x60(%rcx),%xmm10
   0x00000000000020e9 <+8425>:	addps  %xmm13,%xmm12
   0x000000000000211c <+8476>:	addps  %xmm14,%xmm10
   0x0000000000002142 <+8514>:	movaps 0x4a0(%r8),%xmm13
   0x000000000000214a <+8522>:	mulps  0x60(%rcx),%xmm13
   0x000000000000214f <+8527>:	addps  %xmm15,%xmm13
   0x0000000000002153 <+8531>:	movaps 0x4e0(%r8),%xmm15
   0x000000000000215b <+8539>:	mulps  0x60(%rcx),%xmm15
   0x0000000000002182 <+8578>:	addps  %xmm13,%xmm15

697
698	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
699	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
700	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000020a2 <+8354>:	movaps 0x430(%r8),%xmm11
   0x00000000000020c9 <+8393>:	mulps  0x70(%rcx),%xmm11
   0x00000000000020fa <+8442>:	addps  %xmm12,%xmm11
   0x00000000000020fe <+8446>:	movaps 0x470(%r8),%xmm12
   0x0000000000002106 <+8454>:	mulps  0x70(%rcx),%xmm12
   0x000000000000212d <+8493>:	addps  %xmm10,%xmm12
   0x0000000000002131 <+8497>:	movaps 0x4b0(%r8),%xmm10
   0x0000000000002139 <+8505>:	mulps  0x70(%rcx),%xmm10
   0x0000000000002160 <+8544>:	addps  %xmm13,%xmm10
   0x0000000000002175 <+8565>:	movaps 0x4f0(%r8),%xmm14
   0x000000000000217d <+8573>:	mulps  0x70(%rcx),%xmm14
   0x0000000000002196 <+8598>:	addps  %xmm15,%xmm14

701	      }
702
703	      /* A(2,2)*B(2,2) = C(2,2). */
704	      for (i = 0; i < 4; i++)
705	      {
706	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
707	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
708	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002186 <+8582>:	movaps 0x500(%r8),%xmm13
   0x000000000000218e <+8590>:	mulps  0x140(%rcx),%xmm13
   0x00000000000021aa <+8618>:	addps  %xmm11,%xmm13
   0x00000000000021d6 <+8662>:	movaps 0x540(%r8),%xmm11
   0x00000000000021de <+8670>:	mulps  0x140(%rcx),%xmm11
   0x00000000000021fa <+8698>:	addps  %xmm12,%xmm11
   0x0000000000002226 <+8742>:	movaps 0x580(%r8),%xmm12
   0x000000000000222e <+8750>:	mulps  0x140(%rcx),%xmm12
   0x000000000000224a <+8778>:	addps  %xmm10,%xmm12
   0x0000000000002262 <+8802>:	movaps 0x5c0(%r8),%xmm12
   0x000000000000226a <+8810>:	mulps  0x140(%rcx),%xmm12
   0x0000000000002286 <+8838>:	addps  %xmm14,%xmm12

709
710	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
711	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
712	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000021ae <+8622>:	movaps 0x510(%r8),%xmm11
   0x00000000000021b6 <+8630>:	mulps  0x150(%rcx),%xmm11
   0x00000000000021be <+8638>:	addps  %xmm13,%xmm11
   0x00000000000021fe <+8702>:	movaps 0x550(%r8),%xmm12
   0x0000000000002206 <+8710>:	mulps  0x150(%rcx),%xmm12
   0x000000000000220e <+8718>:	addps  %xmm11,%xmm12
   0x000000000000224e <+8782>:	movaps 0x590(%r8),%xmm10
   0x0000000000002256 <+8790>:	mulps  0x150(%rcx),%xmm10
   0x000000000000225e <+8798>:	addps  %xmm12,%xmm10
   0x000000000000228a <+8842>:	movaps 0x5d0(%r8),%xmm14
   0x0000000000002292 <+8850>:	mulps  0x150(%rcx),%xmm14
   0x00000000000022ae <+8878>:	addps  %xmm12,%xmm14

713
714	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
715	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
716	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000219a <+8602>:	movaps 0x520(%r8),%xmm15
   0x00000000000021a2 <+8610>:	mulps  0x160(%rcx),%xmm15
   0x00000000000021d2 <+8658>:	addps  %xmm11,%xmm15
   0x0000000000002212 <+8722>:	movaps 0x560(%r8),%xmm11
   0x000000000000221a <+8730>:	mulps  0x160(%rcx),%xmm11
   0x0000000000002222 <+8738>:	addps  %xmm12,%xmm11
   0x000000000000223a <+8762>:	movaps 0x5a0(%r8),%xmm11
   0x0000000000002242 <+8770>:	mulps  0x160(%rcx),%xmm11
   0x0000000000002272 <+8818>:	addps  %xmm10,%xmm11
   0x00000000000022b2 <+8882>:	movaps 0x5e0(%r8),%xmm12
   0x00000000000022ba <+8890>:	mulps  0x160(%rcx),%xmm12
   0x00000000000022c2 <+8898>:	addps  %xmm14,%xmm12

717
718	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
719	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
720	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000021c2 <+8642>:	movaps 0x530(%r8),%xmm13
   0x00000000000021ca <+8650>:	mulps  0x170(%rcx),%xmm13
   0x00000000000021e6 <+8678>:	addps  %xmm15,%xmm13
   0x00000000000021ea <+8682>:	movaps 0x570(%r8),%xmm15
   0x00000000000021f2 <+8690>:	mulps  0x170(%rcx),%xmm15
   0x0000000000002236 <+8758>:	addps  %xmm11,%xmm15
   0x0000000000002276 <+8822>:	movaps 0x5b0(%r8),%xmm10
   0x000000000000227e <+8830>:	mulps  0x170(%rcx),%xmm10
   0x000000000000229a <+8858>:	addps  %xmm11,%xmm10
   0x000000000000229e <+8862>:	movaps 0x5f0(%r8),%xmm11
   0x00000000000022a6 <+8870>:	mulps  0x170(%rcx),%xmm11
   0x00000000000022d6 <+8918>:	addps  %xmm12,%xmm11

721	      }
722
723	      /* A(2,3)*B(3,2) = C(2,2). */
724	      for (i = 0; i < 4; i++)
725	      {
726	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
727	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
728	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000022c6 <+8902>:	movaps 0x600(%r8),%xmm14
   0x00000000000022ce <+8910>:	mulps  0x240(%rcx),%xmm14
   0x00000000000022ea <+8938>:	addps  %xmm13,%xmm14
   0x0000000000002316 <+8982>:	movaps 0x640(%r8),%xmm13
   0x000000000000231e <+8990>:	mulps  0x240(%rcx),%xmm13
   0x000000000000233a <+9018>:	addps  %xmm15,%xmm13
   0x0000000000002366 <+9062>:	movaps 0x680(%r8),%xmm15
   0x000000000000236e <+9070>:	mulps  0x240(%rcx),%xmm15
   0x000000000000238a <+9098>:	addps  %xmm10,%xmm15
   0x00000000000023b6 <+9142>:	movaps 0x6c0(%r8),%xmm10
   0x00000000000023be <+9150>:	mulps  0x240(%rcx),%xmm10
   0x00000000000023da <+9178>:	addps  %xmm11,%xmm10

729
730	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
731	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
732	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000022ee <+8942>:	movaps 0x610(%r8),%xmm13
   0x00000000000022f6 <+8950>:	mulps  0x250(%rcx),%xmm13
   0x00000000000022fe <+8958>:	addps  %xmm14,%xmm13
   0x000000000000233e <+9022>:	movaps 0x650(%r8),%xmm15
   0x0000000000002346 <+9030>:	mulps  0x250(%rcx),%xmm15
   0x000000000000234e <+9038>:	addps  %xmm13,%xmm15
   0x000000000000238e <+9102>:	movaps 0x690(%r8),%xmm10
   0x0000000000002396 <+9110>:	mulps  0x250(%rcx),%xmm10
   0x000000000000239e <+9118>:	addps  %xmm15,%xmm10
   0x00000000000023de <+9182>:	movaps 0x6d0(%r8),%xmm11
   0x00000000000023e6 <+9190>:	mulps  0x250(%rcx),%xmm11
   0x00000000000023ee <+9198>:	addps  %xmm10,%xmm11

733
734	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
735	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
736	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002302 <+8962>:	movaps 0x620(%r8),%xmm14
   0x000000000000230a <+8970>:	mulps  0x260(%rcx),%xmm14
   0x0000000000002312 <+8978>:	addps  %xmm13,%xmm14
   0x000000000000232a <+9002>:	movaps 0x660(%r8),%xmm14
   0x0000000000002332 <+9010>:	mulps  0x260(%rcx),%xmm14
   0x0000000000002362 <+9058>:	addps  %xmm15,%xmm14
   0x00000000000023a2 <+9122>:	movaps 0x6a0(%r8),%xmm15
   0x00000000000023aa <+9130>:	mulps  0x260(%rcx),%xmm15
   0x00000000000023b2 <+9138>:	addps  %xmm10,%xmm15
   0x00000000000023f2 <+9202>:	movaps 0x6e0(%r8),%xmm10
   0x00000000000023fa <+9210>:	mulps  0x260(%rcx),%xmm10
   0x0000000000002402 <+9218>:	addps  %xmm11,%xmm10

737
738	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
739	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
740	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000022da <+8922>:	movaps 0x630(%r8),%xmm12
   0x00000000000022e2 <+8930>:	mulps  0x270(%rcx),%xmm12
   0x0000000000002326 <+8998>:	addps  %xmm14,%xmm12
   0x0000000000002352 <+9042>:	movaps 0x670(%r8),%xmm13
   0x000000000000235a <+9050>:	mulps  0x270(%rcx),%xmm13
   0x0000000000002376 <+9078>:	addps  %xmm14,%xmm13
   0x000000000000237a <+9082>:	movaps 0x6b0(%r8),%xmm14
   0x0000000000002382 <+9090>:	mulps  0x270(%rcx),%xmm14
   0x00000000000023c6 <+9158>:	addps  %xmm15,%xmm14
   0x00000000000023ca <+9162>:	movaps 0x6f0(%r8),%xmm15
   0x00000000000023d2 <+9170>:	mulps  0x270(%rcx),%xmm15
   0x0000000000002416 <+9238>:	addps  %xmm10,%xmm15

741	      }
742
743	      /* A(2,4)*B(4,2) = C(2,2). */
744	      for (i = 0; i < 4; i++)
745	      {
746	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
747	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
748	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002406 <+9222>:	movaps 0x700(%r8),%xmm11
   0x000000000000240e <+9230>:	mulps  0x340(%rcx),%xmm11
   0x000000000000242a <+9258>:	addps  %xmm12,%xmm11
   0x0000000000002456 <+9302>:	movaps 0x740(%r8),%xmm12
   0x000000000000245e <+9310>:	mulps  0x340(%rcx),%xmm12
   0x000000000000247a <+9338>:	addps  %xmm13,%xmm12
   0x0000000000002492 <+9362>:	movaps 0x780(%r8),%xmm12
   0x000000000000249a <+9370>:	mulps  0x340(%rcx),%xmm12
   0x00000000000024b6 <+9398>:	addps  %xmm14,%xmm12
   0x00000000000024f6 <+9462>:	movaps 0x7c0(%r8),%xmm14
   0x00000000000024fe <+9470>:	mulps  0x340(%rcx),%xmm14
   0x000000000000251a <+9498>:	addps  %xmm15,%xmm14

749
750	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
751	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
752	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000242e <+9262>:	movaps 0x710(%r8),%xmm12
   0x0000000000002436 <+9270>:	mulps  0x350(%rcx),%xmm12
   0x000000000000243e <+9278>:	addps  %xmm11,%xmm12
   0x000000000000247e <+9342>:	movaps 0x750(%r8),%xmm13
   0x0000000000002486 <+9350>:	mulps  0x350(%rcx),%xmm13
   0x000000000000248e <+9358>:	addps  %xmm12,%xmm13
   0x00000000000024ba <+9402>:	movaps 0x790(%r8),%xmm14
   0x00000000000024c2 <+9410>:	mulps  0x350(%rcx),%xmm14
   0x00000000000024de <+9438>:	addps  %xmm12,%xmm14
   0x000000000000251e <+9502>:	movaps 0x7d0(%r8),%xmm15
   0x0000000000002526 <+9510>:	mulps  0x350(%rcx),%xmm15
   0x000000000000252e <+9518>:	addps  %xmm14,%xmm15

753
754	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
755	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
756	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002442 <+9282>:	movaps 0x720(%r8),%xmm11
   0x000000000000244a <+9290>:	mulps  0x360(%rcx),%xmm11
   0x0000000000002452 <+9298>:	addps  %xmm12,%xmm11
   0x000000000000246a <+9322>:	movaps 0x760(%r8),%xmm11
   0x0000000000002472 <+9330>:	mulps  0x360(%rcx),%xmm11
   0x00000000000024a2 <+9378>:	addps  %xmm13,%xmm11
   0x00000000000024ce <+9422>:	movaps 0x7a0(%r8),%xmm11
   0x00000000000024d6 <+9430>:	mulps  0x360(%rcx),%xmm11
   0x00000000000024f2 <+9458>:	addps  %xmm14,%xmm11
   0x0000000000002532 <+9522>:	movaps 0x7e0(%r8),%xmm14
   0x000000000000253a <+9530>:	mulps  0x360(%rcx),%xmm14
   0x0000000000002542 <+9538>:	addps  %xmm15,%xmm14

757
758	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
759	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
760	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000241a <+9242>:	movaps 0x730(%r8),%xmm10
   0x0000000000002422 <+9250>:	mulps  0x370(%rcx),%xmm10
   0x0000000000002466 <+9318>:	addps  %xmm11,%xmm10
   0x00000000000024a6 <+9382>:	movaps 0x770(%r8),%xmm13
   0x00000000000024ae <+9390>:	mulps  0x370(%rcx),%xmm13
   0x00000000000024ca <+9418>:	addps  %xmm11,%xmm13
   0x00000000000024e2 <+9442>:	movaps 0x7b0(%r8),%xmm12
   0x00000000000024ea <+9450>:	mulps  0x370(%rcx),%xmm12
   0x0000000000002506 <+9478>:	addps  %xmm11,%xmm12
   0x000000000000250a <+9482>:	movaps 0x7f0(%r8),%xmm11
   0x0000000000002512 <+9490>:	mulps  0x370(%rcx),%xmm11
   0x0000000000002546 <+9542>:	addps  %xmm14,%xmm11

761	      }
762	    }
763
764	    /* Store C(2,2) block. */
765	    for (i = 0; i < 4; i++)
766	    {
767	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x000000000000254a <+9546>:	mulps  %xmm0,%xmm10
   0x000000000000256c <+9580>:	mulps  %xmm0,%xmm13
   0x0000000000002584 <+9604>:	mulps  %xmm0,%xmm12
   0x000000000000259c <+9628>:	mulps  %xmm0,%xmm11

768	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_22]), C_row[i]);
   0x000000000000254e <+9550>:	addps  0x140(%r9),%xmm10
   0x0000000000002570 <+9584>:	addps  0x150(%r9),%xmm13
   0x0000000000002588 <+9608>:	addps  0x160(%r9),%xmm12
   0x00000000000025a0 <+9632>:	addps  0x170(%r9),%xmm11

769	      _mm_store_ps(&C[i*4+C_OFFSET_22], C_row[i]);
   0x0000000000002560 <+9568>:	movaps %xmm10,0x140(%r9)
   0x0000000000002578 <+9592>:	movaps %xmm13,0x150(%r9)
   0x0000000000002590 <+9616>:	movaps %xmm12,0x160(%r9)
   0x00000000000025a8 <+9640>:	movaps %xmm11,0x170(%r9)

770	    }
771
772	    /* Reset C(2,3) matrix accumulators */
773	    C_row[0] = _mm_setzero_ps();
   0x0000000000002568 <+9576>:	xorps  %xmm10,%xmm10

774	    C_row[1] = _mm_setzero_ps();
   0x0000000000002580 <+9600>:	xorps  %xmm13,%xmm13

775	    C_row[2] = _mm_setzero_ps();
   0x0000000000002598 <+9624>:	xorps  %xmm12,%xmm12

776	    C_row[3] = _mm_setzero_ps();
   0x00000000000025b0 <+9648>:	xorps  %xmm11,%xmm11

777
778	    if (norm_product[4][2] &&
   0x0000000000002556 <+9558>:	movss  0x148(%rsp),%xmm14
   0x00000000000025b4 <+9652>:	comiss %xmm1,%xmm14
   0x00000000000025b8 <+9656>:	jb     0x2aea <stream_kernel+10986>

779	        norm_product[5][6] &&
   0x00000000000025be <+9662>:	movss  0xe8(%rsp),%xmm14
   0x00000000000025c8 <+9672>:	comiss %xmm1,%xmm14
   0x00000000000025cc <+9676>:	jb     0x2aea <stream_kernel+10986>

780	        norm_product[6][10] &&
   0x00000000000025d2 <+9682>:	movss  0xe0(%rsp),%xmm14
   0x00000000000025dc <+9692>:	comiss %xmm1,%xmm14
   0x00000000000025e0 <+9696>:	jb     0x2aea <stream_kernel+10986>

781	        norm_product[7][14])
   0x00000000000025e6 <+9702>:	movss  0xd8(%rsp),%xmm14
   0x00000000000025f0 <+9712>:	comiss %xmm1,%xmm14
   0x00000000000025f4 <+9716>:	jb     0x2aea <stream_kernel+10986>

782	    {
783	      /* A(2,1)*B(1,3) = C(2,3). */
784	      for (i = 0; i < 4; i++)
785	      {
786	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
787	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
788	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000025fa <+9722>:	movaps 0x400(%r8),%xmm10
   0x000000000000261a <+9754>:	movaps 0x440(%r8),%xmm15
   0x000000000000262a <+9770>:	mulps  0x80(%rcx),%xmm10
   0x000000000000264a <+9802>:	mulps  0x80(%rcx),%xmm15
   0x0000000000002672 <+9842>:	movaps 0x480(%r8),%xmm13
   0x000000000000267a <+9850>:	mulps  0x80(%rcx),%xmm13
   0x00000000000026ae <+9902>:	movaps 0x4c0(%r8),%xmm14
   0x00000000000026b6 <+9910>:	mulps  0x80(%rcx),%xmm14

789
790	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
791	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
792	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002602 <+9730>:	movaps 0x410(%r8),%xmm13
   0x0000000000002622 <+9762>:	movaps 0x450(%r8),%xmm14
   0x0000000000002632 <+9778>:	mulps  0x90(%rcx),%xmm13
   0x0000000000002652 <+9810>:	mulps  0x90(%rcx),%xmm14
   0x000000000000265a <+9818>:	addps  %xmm10,%xmm13
   0x0000000000002696 <+9878>:	addps  %xmm15,%xmm14
   0x000000000000269a <+9882>:	movaps 0x490(%r8),%xmm15
   0x00000000000026a2 <+9890>:	mulps  0x90(%rcx),%xmm15
   0x00000000000026d2 <+9938>:	addps  %xmm13,%xmm15
   0x00000000000026fe <+9982>:	movaps 0x4d0(%r8),%xmm13
   0x0000000000002706 <+9990>:	mulps  0x90(%rcx),%xmm13
   0x000000000000270e <+9998>:	addps  %xmm14,%xmm13

793
794	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
795	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
796	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000260a <+9738>:	movaps 0x420(%r8),%xmm12
   0x000000000000263a <+9786>:	mulps  0xa0(%rcx),%xmm12
   0x000000000000265e <+9822>:	movaps 0x460(%r8),%xmm10
   0x0000000000002666 <+9830>:	mulps  0xa0(%rcx),%xmm10
   0x000000000000266e <+9838>:	addps  %xmm13,%xmm12
   0x00000000000026aa <+9898>:	addps  %xmm14,%xmm10
   0x00000000000026d6 <+9942>:	movaps 0x4a0(%r8),%xmm13
   0x00000000000026de <+9950>:	mulps  0xa0(%rcx),%xmm13
   0x00000000000026e6 <+9958>:	addps  %xmm15,%xmm13
   0x00000000000026ea <+9962>:	movaps 0x4e0(%r8),%xmm15
   0x00000000000026f2 <+9970>:	mulps  0xa0(%rcx),%xmm15
   0x0000000000002722 <+10018>:	addps  %xmm13,%xmm15

797
798	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
799	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
800	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002612 <+9746>:	movaps 0x430(%r8),%xmm11
   0x0000000000002642 <+9794>:	mulps  0xb0(%rcx),%xmm11
   0x0000000000002682 <+9858>:	addps  %xmm12,%xmm11
   0x0000000000002686 <+9862>:	movaps 0x470(%r8),%xmm12
   0x000000000000268e <+9870>:	mulps  0xb0(%rcx),%xmm12
   0x00000000000026be <+9918>:	addps  %xmm10,%xmm12
   0x00000000000026c2 <+9922>:	movaps 0x4b0(%r8),%xmm10
   0x00000000000026ca <+9930>:	mulps  0xb0(%rcx),%xmm10
   0x00000000000026fa <+9978>:	addps  %xmm13,%xmm10
   0x0000000000002712 <+10002>:	movaps 0x4f0(%r8),%xmm14
   0x000000000000271a <+10010>:	mulps  0xb0(%rcx),%xmm14
   0x0000000000002736 <+10038>:	addps  %xmm15,%xmm14

801	      }
802
803	      /* A(2,2)*B(2,3) = C(2,3). */
804	      for (i = 0; i < 4; i++)
805	      {
806	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
807	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
808	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002726 <+10022>:	movaps 0x500(%r8),%xmm13
   0x000000000000272e <+10030>:	mulps  0x180(%rcx),%xmm13
   0x000000000000274a <+10058>:	addps  %xmm11,%xmm13
   0x0000000000002776 <+10102>:	movaps 0x540(%r8),%xmm11
   0x000000000000277e <+10110>:	mulps  0x180(%rcx),%xmm11
   0x000000000000279a <+10138>:	addps  %xmm12,%xmm11
   0x00000000000027c6 <+10182>:	movaps 0x580(%r8),%xmm12
   0x00000000000027ce <+10190>:	mulps  0x180(%rcx),%xmm12
   0x00000000000027ea <+10218>:	addps  %xmm10,%xmm12
   0x0000000000002802 <+10242>:	movaps 0x5c0(%r8),%xmm12
   0x000000000000280a <+10250>:	mulps  0x180(%rcx),%xmm12
   0x0000000000002826 <+10278>:	addps  %xmm14,%xmm12

809
810	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
811	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
812	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000274e <+10062>:	movaps 0x510(%r8),%xmm11
   0x0000000000002756 <+10070>:	mulps  0x190(%rcx),%xmm11
   0x000000000000275e <+10078>:	addps  %xmm13,%xmm11
   0x000000000000279e <+10142>:	movaps 0x550(%r8),%xmm12
   0x00000000000027a6 <+10150>:	mulps  0x190(%rcx),%xmm12
   0x00000000000027ae <+10158>:	addps  %xmm11,%xmm12
   0x00000000000027ee <+10222>:	movaps 0x590(%r8),%xmm10
   0x00000000000027f6 <+10230>:	mulps  0x190(%rcx),%xmm10
   0x00000000000027fe <+10238>:	addps  %xmm12,%xmm10
   0x000000000000282a <+10282>:	movaps 0x5d0(%r8),%xmm14
   0x0000000000002832 <+10290>:	mulps  0x190(%rcx),%xmm14
   0x000000000000284e <+10318>:	addps  %xmm12,%xmm14

813
814	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
815	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
816	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000273a <+10042>:	movaps 0x520(%r8),%xmm15
   0x0000000000002742 <+10050>:	mulps  0x1a0(%rcx),%xmm15
   0x0000000000002772 <+10098>:	addps  %xmm11,%xmm15
   0x00000000000027b2 <+10162>:	movaps 0x560(%r8),%xmm11
   0x00000000000027ba <+10170>:	mulps  0x1a0(%rcx),%xmm11
   0x00000000000027c2 <+10178>:	addps  %xmm12,%xmm11
   0x00000000000027da <+10202>:	movaps 0x5a0(%r8),%xmm11
   0x00000000000027e2 <+10210>:	mulps  0x1a0(%rcx),%xmm11
   0x0000000000002812 <+10258>:	addps  %xmm10,%xmm11
   0x0000000000002852 <+10322>:	movaps 0x5e0(%r8),%xmm12
   0x000000000000285a <+10330>:	mulps  0x1a0(%rcx),%xmm12
   0x0000000000002862 <+10338>:	addps  %xmm14,%xmm12

817
818	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
819	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
820	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002762 <+10082>:	movaps 0x530(%r8),%xmm13
   0x000000000000276a <+10090>:	mulps  0x1b0(%rcx),%xmm13
   0x0000000000002786 <+10118>:	addps  %xmm15,%xmm13
   0x000000000000278a <+10122>:	movaps 0x570(%r8),%xmm15
   0x0000000000002792 <+10130>:	mulps  0x1b0(%rcx),%xmm15
   0x00000000000027d6 <+10198>:	addps  %xmm11,%xmm15
   0x0000000000002816 <+10262>:	movaps 0x5b0(%r8),%xmm10
   0x000000000000281e <+10270>:	mulps  0x1b0(%rcx),%xmm10
   0x000000000000283a <+10298>:	addps  %xmm11,%xmm10
   0x000000000000283e <+10302>:	movaps 0x5f0(%r8),%xmm11
   0x0000000000002846 <+10310>:	mulps  0x1b0(%rcx),%xmm11
   0x0000000000002876 <+10358>:	addps  %xmm12,%xmm11

821	      }
822
823	      /* A(2,3)*B(3,3) = C(2,3). */
824	      for (i = 0; i < 4; i++)
825	      {
826	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
827	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
828	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002866 <+10342>:	movaps 0x600(%r8),%xmm14
   0x000000000000286e <+10350>:	mulps  0x280(%rcx),%xmm14
   0x000000000000288a <+10378>:	addps  %xmm13,%xmm14
   0x00000000000028b6 <+10422>:	movaps 0x640(%r8),%xmm13
   0x00000000000028be <+10430>:	mulps  0x280(%rcx),%xmm13
   0x00000000000028da <+10458>:	addps  %xmm15,%xmm13
   0x0000000000002906 <+10502>:	movaps 0x680(%r8),%xmm15
   0x000000000000290e <+10510>:	mulps  0x280(%rcx),%xmm15
   0x000000000000292a <+10538>:	addps  %xmm10,%xmm15
   0x0000000000002956 <+10582>:	movaps 0x6c0(%r8),%xmm10
   0x000000000000295e <+10590>:	mulps  0x280(%rcx),%xmm10
   0x000000000000297a <+10618>:	addps  %xmm11,%xmm10

829
830	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
831	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
832	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000288e <+10382>:	movaps 0x610(%r8),%xmm13
   0x0000000000002896 <+10390>:	mulps  0x290(%rcx),%xmm13
   0x000000000000289e <+10398>:	addps  %xmm14,%xmm13
   0x00000000000028de <+10462>:	movaps 0x650(%r8),%xmm15
   0x00000000000028e6 <+10470>:	mulps  0x290(%rcx),%xmm15
   0x00000000000028ee <+10478>:	addps  %xmm13,%xmm15
   0x000000000000292e <+10542>:	movaps 0x690(%r8),%xmm10
   0x0000000000002936 <+10550>:	mulps  0x290(%rcx),%xmm10
   0x000000000000293e <+10558>:	addps  %xmm15,%xmm10
   0x000000000000297e <+10622>:	movaps 0x6d0(%r8),%xmm11
   0x0000000000002986 <+10630>:	mulps  0x290(%rcx),%xmm11
   0x000000000000298e <+10638>:	addps  %xmm10,%xmm11

833
834	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
835	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
836	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000028a2 <+10402>:	movaps 0x620(%r8),%xmm14
   0x00000000000028aa <+10410>:	mulps  0x2a0(%rcx),%xmm14
   0x00000000000028b2 <+10418>:	addps  %xmm13,%xmm14
   0x00000000000028ca <+10442>:	movaps 0x660(%r8),%xmm14
   0x00000000000028d2 <+10450>:	mulps  0x2a0(%rcx),%xmm14
   0x0000000000002902 <+10498>:	addps  %xmm15,%xmm14
   0x0000000000002942 <+10562>:	movaps 0x6a0(%r8),%xmm15
   0x000000000000294a <+10570>:	mulps  0x2a0(%rcx),%xmm15
   0x0000000000002952 <+10578>:	addps  %xmm10,%xmm15
   0x0000000000002992 <+10642>:	movaps 0x6e0(%r8),%xmm10
   0x000000000000299a <+10650>:	mulps  0x2a0(%rcx),%xmm10
   0x00000000000029a2 <+10658>:	addps  %xmm11,%xmm10

837
838	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
839	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
840	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000287a <+10362>:	movaps 0x630(%r8),%xmm12
   0x0000000000002882 <+10370>:	mulps  0x2b0(%rcx),%xmm12
   0x00000000000028c6 <+10438>:	addps  %xmm14,%xmm12
   0x00000000000028f2 <+10482>:	movaps 0x670(%r8),%xmm13
   0x00000000000028fa <+10490>:	mulps  0x2b0(%rcx),%xmm13
   0x0000000000002916 <+10518>:	addps  %xmm14,%xmm13
   0x000000000000291a <+10522>:	movaps 0x6b0(%r8),%xmm14
   0x0000000000002922 <+10530>:	mulps  0x2b0(%rcx),%xmm14
   0x0000000000002966 <+10598>:	addps  %xmm15,%xmm14
   0x000000000000296a <+10602>:	movaps 0x6f0(%r8),%xmm15
   0x0000000000002972 <+10610>:	mulps  0x2b0(%rcx),%xmm15
   0x00000000000029b6 <+10678>:	addps  %xmm10,%xmm15

841	      }
842
843	      /* A(2,4)*B(4,3) = C(2,3). */
844	      for (i = 0; i < 4; i++)
845	      {
846	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
847	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
848	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000029a6 <+10662>:	movaps 0x700(%r8),%xmm11
   0x00000000000029ae <+10670>:	mulps  0x380(%rcx),%xmm11
   0x00000000000029ca <+10698>:	addps  %xmm12,%xmm11
   0x00000000000029f6 <+10742>:	movaps 0x740(%r8),%xmm12
   0x00000000000029fe <+10750>:	mulps  0x380(%rcx),%xmm12
   0x0000000000002a1a <+10778>:	addps  %xmm13,%xmm12
   0x0000000000002a32 <+10802>:	movaps 0x780(%r8),%xmm12
   0x0000000000002a3a <+10810>:	mulps  0x380(%rcx),%xmm12
   0x0000000000002a56 <+10838>:	addps  %xmm14,%xmm12
   0x0000000000002a96 <+10902>:	movaps 0x7c0(%r8),%xmm14
   0x0000000000002a9e <+10910>:	mulps  0x380(%rcx),%xmm14
   0x0000000000002aba <+10938>:	addps  %xmm15,%xmm14

849
850	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
851	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
852	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000029ce <+10702>:	movaps 0x710(%r8),%xmm12
   0x00000000000029d6 <+10710>:	mulps  0x390(%rcx),%xmm12
   0x00000000000029de <+10718>:	addps  %xmm11,%xmm12
   0x0000000000002a1e <+10782>:	movaps 0x750(%r8),%xmm13
   0x0000000000002a26 <+10790>:	mulps  0x390(%rcx),%xmm13
   0x0000000000002a2e <+10798>:	addps  %xmm12,%xmm13
   0x0000000000002a5a <+10842>:	movaps 0x790(%r8),%xmm14
   0x0000000000002a62 <+10850>:	mulps  0x390(%rcx),%xmm14
   0x0000000000002a7e <+10878>:	addps  %xmm12,%xmm14
   0x0000000000002abe <+10942>:	movaps 0x7d0(%r8),%xmm15
   0x0000000000002ac6 <+10950>:	mulps  0x390(%rcx),%xmm15
   0x0000000000002ace <+10958>:	addps  %xmm14,%xmm15

853
854	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
855	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
856	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000029e2 <+10722>:	movaps 0x720(%r8),%xmm11
   0x00000000000029ea <+10730>:	mulps  0x3a0(%rcx),%xmm11
   0x00000000000029f2 <+10738>:	addps  %xmm12,%xmm11
   0x0000000000002a0a <+10762>:	movaps 0x760(%r8),%xmm11
   0x0000000000002a12 <+10770>:	mulps  0x3a0(%rcx),%xmm11
   0x0000000000002a42 <+10818>:	addps  %xmm13,%xmm11
   0x0000000000002a6e <+10862>:	movaps 0x7a0(%r8),%xmm11
   0x0000000000002a76 <+10870>:	mulps  0x3a0(%rcx),%xmm11
   0x0000000000002a92 <+10898>:	addps  %xmm14,%xmm11
   0x0000000000002ad2 <+10962>:	movaps 0x7e0(%r8),%xmm14
   0x0000000000002ada <+10970>:	mulps  0x3a0(%rcx),%xmm14
   0x0000000000002ae2 <+10978>:	addps  %xmm15,%xmm14

857
858	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
859	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
860	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000029ba <+10682>:	movaps 0x730(%r8),%xmm10
   0x00000000000029c2 <+10690>:	mulps  0x3b0(%rcx),%xmm10
   0x0000000000002a06 <+10758>:	addps  %xmm11,%xmm10
   0x0000000000002a46 <+10822>:	movaps 0x770(%r8),%xmm13
   0x0000000000002a4e <+10830>:	mulps  0x3b0(%rcx),%xmm13
   0x0000000000002a6a <+10858>:	addps  %xmm11,%xmm13
   0x0000000000002a82 <+10882>:	movaps 0x7b0(%r8),%xmm12
   0x0000000000002a8a <+10890>:	mulps  0x3b0(%rcx),%xmm12
   0x0000000000002aa6 <+10918>:	addps  %xmm11,%xmm12
   0x0000000000002aaa <+10922>:	movaps 0x7f0(%r8),%xmm11
   0x0000000000002ab2 <+10930>:	mulps  0x3b0(%rcx),%xmm11
   0x0000000000002ae6 <+10982>:	addps  %xmm14,%xmm11

861	      }
862	    }
863
864	    /* Store C(2,3) block. */
865	    for (i = 0; i < 4; i++)
866	    {
867	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000002aea <+10986>:	mulps  %xmm0,%xmm10
   0x0000000000002b02 <+11010>:	mulps  %xmm0,%xmm13
   0x0000000000002b1a <+11034>:	mulps  %xmm0,%xmm12
   0x0000000000002b32 <+11058>:	mulps  %xmm0,%xmm11

868	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_23]), C_row[i]);
   0x0000000000002aee <+10990>:	addps  0x180(%r9),%xmm10
   0x0000000000002b06 <+11014>:	addps  0x190(%r9),%xmm13
   0x0000000000002b1e <+11038>:	addps  0x1a0(%r9),%xmm12
   0x0000000000002b36 <+11062>:	addps  0x1b0(%r9),%xmm11

869	      _mm_store_ps(&C[i*4+C_OFFSET_23], C_row[i]);
   0x0000000000002af6 <+10998>:	movaps %xmm10,0x180(%r9)
   0x0000000000002b0e <+11022>:	movaps %xmm13,0x190(%r9)
   0x0000000000002b26 <+11046>:	movaps %xmm12,0x1a0(%r9)
   0x0000000000002b3e <+11070>:	movaps %xmm11,0x1b0(%r9)

870	    }
871
872	    /* Reset C(2,4) matrix accumulators */
873	    C_row[0] = _mm_setzero_ps();
   0x0000000000002afe <+11006>:	xorps  %xmm10,%xmm10

874	    C_row[1] = _mm_setzero_ps();
   0x0000000000002b16 <+11030>:	xorps  %xmm13,%xmm13

875	    C_row[2] = _mm_setzero_ps();
   0x0000000000002b2e <+11054>:	xorps  %xmm12,%xmm12

876	    C_row[3] = _mm_setzero_ps();
   0x0000000000002b46 <+11078>:	xorps  %xmm11,%xmm11

877
878	    if (norm_product[4][3] &&
   0x0000000000002b4a <+11082>:	comiss %xmm1,%xmm9
   0x0000000000002b4e <+11086>:	jb     0x3080 <stream_kernel+12416>

879	        norm_product[5][7] &&
   0x0000000000002b54 <+11092>:	movss  0xd0(%rsp),%xmm9
   0x0000000000002b5e <+11102>:	comiss %xmm1,%xmm9
   0x0000000000002b62 <+11106>:	jb     0x3080 <stream_kernel+12416>

880	        norm_product[6][11] &&
   0x0000000000002b68 <+11112>:	movss  0xc8(%rsp),%xmm9
   0x0000000000002b72 <+11122>:	comiss %xmm1,%xmm9
   0x0000000000002b76 <+11126>:	jb     0x3080 <stream_kernel+12416>

881	        norm_product[7][15])
   0x0000000000002b7c <+11132>:	movss  0xb8(%rsp),%xmm9
   0x0000000000002b86 <+11142>:	comiss %xmm1,%xmm9
   0x0000000000002b8a <+11146>:	jb     0x3080 <stream_kernel+12416>

882	    {
883	      /* A(2,1)*B(1,4) = C(2,4). */
884	      for (i = 0; i < 4; i++)
885	      {
886	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
887	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
888	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b90 <+11152>:	movaps 0x400(%r8),%xmm14
   0x0000000000002bb0 <+11184>:	movaps 0x440(%r8),%xmm11
   0x0000000000002bc8 <+11208>:	mulps  0xc0(%rcx),%xmm14
   0x0000000000002be8 <+11240>:	mulps  0xc0(%rcx),%xmm11
   0x0000000000002c04 <+11268>:	movaps 0x4c0(%r8),%xmm14
   0x0000000000002c0c <+11276>:	mulps  0xc0(%rcx),%xmm14
   0x0000000000002c18 <+11288>:	movaps 0x480(%r8),%xmm10
   0x0000000000002c20 <+11296>:	mulps  0xc0(%rcx),%xmm10

889
890	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
891	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
892	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b98 <+11160>:	movaps 0x410(%r8),%xmm10
   0x0000000000002bb8 <+11192>:	movaps 0x450(%r8),%xmm9
   0x0000000000002bd0 <+11216>:	mulps  0xd0(%rcx),%xmm10
   0x0000000000002bf0 <+11248>:	mulps  0xd0(%rcx),%xmm9
   0x0000000000002c00 <+11264>:	addps  %xmm14,%xmm10
   0x0000000000002c3c <+11324>:	addps  %xmm11,%xmm9
   0x0000000000002c40 <+11328>:	movaps 0x490(%r8),%xmm11
   0x0000000000002c48 <+11336>:	mulps  0xd0(%rcx),%xmm11
   0x0000000000002c78 <+11384>:	addps  %xmm10,%xmm11
   0x0000000000002c7c <+11388>:	movaps 0x4d0(%r8),%xmm10
   0x0000000000002c84 <+11396>:	mulps  0xd0(%rcx),%xmm10
   0x0000000000002ca0 <+11424>:	addps  %xmm14,%xmm10

893
894	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
895	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
896	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002ba0 <+11168>:	movaps 0x420(%r8),%xmm12
   0x0000000000002bc0 <+11200>:	movaps 0x460(%r8),%xmm15
   0x0000000000002bd8 <+11224>:	mulps  0xe0(%rcx),%xmm12
   0x0000000000002bf8 <+11256>:	mulps  0xe0(%rcx),%xmm15
   0x0000000000002c14 <+11284>:	addps  %xmm10,%xmm12
   0x0000000000002c50 <+11344>:	addps  %xmm9,%xmm15
   0x0000000000002c54 <+11348>:	movaps 0x4a0(%r8),%xmm9
   0x0000000000002c5c <+11356>:	mulps  0xe0(%rcx),%xmm9
   0x0000000000002c68 <+11368>:	movaps 0x4e0(%r8),%xmm15
   0x0000000000002c70 <+11376>:	mulps  0xe0(%rcx),%xmm15
   0x0000000000002c8c <+11404>:	addps  %xmm11,%xmm9
   0x0000000000002cc8 <+11464>:	addps  %xmm10,%xmm15

897
898	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
899	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
900	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002ba8 <+11176>:	movaps 0x430(%r8),%xmm13
   0x0000000000002be0 <+11232>:	mulps  0xf0(%rcx),%xmm13
   0x0000000000002c28 <+11304>:	addps  %xmm12,%xmm13
   0x0000000000002c2c <+11308>:	movaps 0x470(%r8),%xmm12
   0x0000000000002c34 <+11316>:	mulps  0xf0(%rcx),%xmm12
   0x0000000000002c64 <+11364>:	addps  %xmm15,%xmm12
   0x0000000000002c90 <+11408>:	movaps 0x4b0(%r8),%xmm11
   0x0000000000002c98 <+11416>:	mulps  0xf0(%rcx),%xmm11
   0x0000000000002cb4 <+11444>:	addps  %xmm9,%xmm11
   0x0000000000002cb8 <+11448>:	movaps 0x4f0(%r8),%xmm9
   0x0000000000002cc0 <+11456>:	mulps  0xf0(%rcx),%xmm9
   0x0000000000002cdc <+11484>:	addps  %xmm15,%xmm9

901	      }
902
903	      /* A(2,2)*B(2,4) = C(2,4). */
904	      for (i = 0; i < 4; i++)
905	      {
906	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
907	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
908	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002ccc <+11468>:	movaps 0x500(%r8),%xmm10
   0x0000000000002cd4 <+11476>:	mulps  0x1c0(%rcx),%xmm10
   0x0000000000002ce0 <+11488>:	movaps 0x540(%r8),%xmm15
   0x0000000000002ce8 <+11496>:	mulps  0x1c0(%rcx),%xmm15
   0x0000000000002cf0 <+11504>:	addps  %xmm13,%xmm10
   0x0000000000002d04 <+11524>:	addps  %xmm12,%xmm15
   0x0000000000002d44 <+11588>:	movaps 0x580(%r8),%xmm14
   0x0000000000002d4c <+11596>:	mulps  0x1c0(%rcx),%xmm14
   0x0000000000002d6c <+11628>:	movaps 0x5c0(%r8),%xmm12
   0x0000000000002d74 <+11636>:	mulps  0x1c0(%rcx),%xmm12
   0x0000000000002d90 <+11664>:	addps  %xmm11,%xmm14
   0x0000000000002da4 <+11684>:	addps  %xmm9,%xmm12

909
910	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
911	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
912	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002cf4 <+11508>:	movaps 0x510(%r8),%xmm13
   0x0000000000002cfc <+11516>:	mulps  0x1d0(%rcx),%xmm13
   0x0000000000002d08 <+11528>:	movaps 0x550(%r8),%xmm12
   0x0000000000002d10 <+11536>:	mulps  0x1d0(%rcx),%xmm12
   0x0000000000002d18 <+11544>:	addps  %xmm10,%xmm13
   0x0000000000002d54 <+11604>:	addps  %xmm15,%xmm12
   0x0000000000002d94 <+11668>:	movaps 0x590(%r8),%xmm11
   0x0000000000002d9c <+11676>:	mulps  0x1d0(%rcx),%xmm11
   0x0000000000002da8 <+11688>:	movaps 0x5d0(%r8),%xmm9
   0x0000000000002db0 <+11696>:	mulps  0x1d0(%rcx),%xmm9
   0x0000000000002db8 <+11704>:	addps  %xmm14,%xmm11
   0x0000000000002df4 <+11764>:	addps  %xmm12,%xmm9

913
914	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
915	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
916	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002ca4 <+11428>:	movaps 0x520(%r8),%xmm14
   0x0000000000002cac <+11436>:	mulps  0x1e0(%rcx),%xmm14
   0x0000000000002d2c <+11564>:	addps  %xmm13,%xmm14
   0x0000000000002d30 <+11568>:	movaps 0x560(%r8),%xmm13
   0x0000000000002d38 <+11576>:	mulps  0x1e0(%rcx),%xmm13
   0x0000000000002d68 <+11624>:	addps  %xmm12,%xmm13
   0x0000000000002d80 <+11648>:	movaps 0x5a0(%r8),%xmm13
   0x0000000000002d88 <+11656>:	mulps  0x1e0(%rcx),%xmm13
   0x0000000000002dcc <+11724>:	addps  %xmm11,%xmm13
   0x0000000000002df8 <+11768>:	movaps 0x5e0(%r8),%xmm12
   0x0000000000002e00 <+11776>:	mulps  0x1e0(%rcx),%xmm12
   0x0000000000002e1c <+11804>:	addps  %xmm9,%xmm12

917
918	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
919	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
920	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002d1c <+11548>:	movaps 0x530(%r8),%xmm10
   0x0000000000002d24 <+11556>:	mulps  0x1f0(%rcx),%xmm10
   0x0000000000002d40 <+11584>:	addps  %xmm14,%xmm10
   0x0000000000002d58 <+11608>:	movaps 0x570(%r8),%xmm15
   0x0000000000002d60 <+11616>:	mulps  0x1f0(%rcx),%xmm15
   0x0000000000002d7c <+11644>:	addps  %xmm13,%xmm15
   0x0000000000002dbc <+11708>:	movaps 0x5b0(%r8),%xmm14
   0x0000000000002dc4 <+11716>:	mulps  0x1f0(%rcx),%xmm14
   0x0000000000002de0 <+11744>:	addps  %xmm13,%xmm14
   0x0000000000002de4 <+11748>:	movaps 0x5f0(%r8),%xmm13
   0x0000000000002dec <+11756>:	mulps  0x1f0(%rcx),%xmm13
   0x0000000000002e30 <+11824>:	addps  %xmm12,%xmm13

921	      }
922
923	      /* A(2,3)*B(3,4) = C(2,4). */
924	      for (i = 0; i < 4; i++)
925	      {
926	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
927	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
928	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002dd0 <+11728>:	movaps 0x600(%r8),%xmm11
   0x0000000000002dd8 <+11736>:	mulps  0x2c0(%rcx),%xmm11
   0x0000000000002e08 <+11784>:	addps  %xmm10,%xmm11
   0x0000000000002e20 <+11808>:	movaps 0x640(%r8),%xmm9
   0x0000000000002e28 <+11816>:	mulps  0x2c0(%rcx),%xmm9
   0x0000000000002e6c <+11884>:	addps  %xmm15,%xmm9
   0x0000000000002e98 <+11928>:	movaps 0x680(%r8),%xmm9
   0x0000000000002ea0 <+11936>:	mulps  0x2c0(%rcx),%xmm9
   0x0000000000002eac <+11948>:	movaps 0x6c0(%r8),%xmm15
   0x0000000000002eb4 <+11956>:	mulps  0x2c0(%rcx),%xmm15
   0x0000000000002ed0 <+11984>:	addps  %xmm14,%xmm9
   0x0000000000002ee4 <+12004>:	addps  %xmm13,%xmm15

929
930	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
931	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
932	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002e0c <+11788>:	movaps 0x610(%r8),%xmm10
   0x0000000000002e14 <+11796>:	mulps  0x2d0(%rcx),%xmm10
   0x0000000000002e44 <+11844>:	addps  %xmm11,%xmm10
   0x0000000000002e70 <+11888>:	movaps 0x650(%r8),%xmm15
   0x0000000000002e78 <+11896>:	mulps  0x2d0(%rcx),%xmm15
   0x0000000000002e94 <+11924>:	addps  %xmm9,%xmm15
   0x0000000000002ed4 <+11988>:	movaps 0x690(%r8),%xmm14
   0x0000000000002edc <+11996>:	mulps  0x2d0(%rcx),%xmm14
   0x0000000000002ee8 <+12008>:	movaps 0x6d0(%r8),%xmm13
   0x0000000000002ef0 <+12016>:	mulps  0x2d0(%rcx),%xmm13
   0x0000000000002ef8 <+12024>:	addps  %xmm9,%xmm14
   0x0000000000002f0c <+12044>:	addps  %xmm15,%xmm13

933
934	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
935	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
936	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002e34 <+11828>:	movaps 0x620(%r8),%xmm12
   0x0000000000002e3c <+11836>:	mulps  0x2e0(%rcx),%xmm12
   0x0000000000002e48 <+11848>:	movaps 0x660(%r8),%xmm11
   0x0000000000002e50 <+11856>:	mulps  0x2e0(%rcx),%xmm11
   0x0000000000002e58 <+11864>:	addps  %xmm10,%xmm12
   0x0000000000002ea8 <+11944>:	addps  %xmm15,%xmm11
   0x0000000000002efc <+12028>:	movaps 0x6a0(%r8),%xmm9
   0x0000000000002f04 <+12036>:	mulps  0x2e0(%rcx),%xmm9
   0x0000000000002f20 <+12064>:	addps  %xmm14,%xmm9
   0x0000000000002f24 <+12068>:	movaps 0x6e0(%r8),%xmm14
   0x0000000000002f2c <+12076>:	mulps  0x2e0(%rcx),%xmm14
   0x0000000000002f48 <+12104>:	addps  %xmm13,%xmm14

937
938	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
939	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
940	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002e5c <+11868>:	movaps 0x630(%r8),%xmm10
   0x0000000000002e64 <+11876>:	mulps  0x2f0(%rcx),%xmm10
   0x0000000000002e80 <+11904>:	addps  %xmm12,%xmm10
   0x0000000000002e84 <+11908>:	movaps 0x670(%r8),%xmm12
   0x0000000000002e8c <+11916>:	mulps  0x2f0(%rcx),%xmm12
   0x0000000000002ebc <+11964>:	addps  %xmm11,%xmm12
   0x0000000000002ec0 <+11968>:	movaps 0x6b0(%r8),%xmm11
   0x0000000000002ec8 <+11976>:	mulps  0x2f0(%rcx),%xmm11
   0x0000000000002f34 <+12084>:	addps  %xmm9,%xmm11
   0x0000000000002f38 <+12088>:	movaps 0x6f0(%r8),%xmm9
   0x0000000000002f40 <+12096>:	mulps  0x2f0(%rcx),%xmm9
   0x0000000000002f5c <+12124>:	addps  %xmm14,%xmm9

941	      }
942
943	      /* A(2,4)*B(4,4) = C(2,4). */
944	      for (i = 0; i < 4; i++)
945	      {
946	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
947	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
948	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002f10 <+12048>:	movaps 0x700(%r8),%xmm15
   0x0000000000002f18 <+12056>:	mulps  0x3c0(%rcx),%xmm15
   0x0000000000002f60 <+12128>:	movaps 0x740(%r8),%xmm14
   0x0000000000002f68 <+12136>:	mulps  0x3c0(%rcx),%xmm14
   0x0000000000002f70 <+12144>:	addps  %xmm10,%xmm15
   0x0000000000002f84 <+12164>:	addps  %xmm12,%xmm14
   0x0000000000002f9c <+12188>:	movaps 0x780(%r8),%xmm15
   0x0000000000002fa4 <+12196>:	mulps  0x3c0(%rcx),%xmm15
   0x0000000000003010 <+12304>:	addps  %xmm11,%xmm15
   0x0000000000003028 <+12328>:	movaps 0x7c0(%r8),%xmm15
   0x0000000000003030 <+12336>:	mulps  0x3c0(%rcx),%xmm15
   0x0000000000003060 <+12384>:	addps  %xmm9,%xmm15

949
950	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
951	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
952	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002f74 <+12148>:	movaps 0x710(%r8),%xmm10
   0x0000000000002f7c <+12156>:	mulps  0x3d0(%rcx),%xmm10
   0x0000000000002f88 <+12168>:	movaps 0x750(%r8),%xmm12
   0x0000000000002f90 <+12176>:	mulps  0x3d0(%rcx),%xmm12
   0x0000000000002f98 <+12184>:	addps  %xmm15,%xmm10
   0x0000000000002fc0 <+12224>:	addps  %xmm14,%xmm12
   0x0000000000003014 <+12308>:	movaps 0x790(%r8),%xmm11
   0x000000000000301c <+12316>:	mulps  0x3d0(%rcx),%xmm11
   0x0000000000003024 <+12324>:	addps  %xmm15,%xmm11
   0x0000000000003064 <+12388>:	movaps 0x7d0(%r8),%xmm9
   0x000000000000306c <+12396>:	mulps  0x3d0(%rcx),%xmm9
   0x0000000000003074 <+12404>:	addps  %xmm15,%xmm9

953
954	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
955	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
956	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002f4c <+12108>:	movaps 0x720(%r8),%xmm13
   0x0000000000002f54 <+12116>:	mulps  0x3e0(%rcx),%xmm13
   0x0000000000002fac <+12204>:	addps  %xmm10,%xmm13
   0x0000000000002fc4 <+12228>:	movaps 0x760(%r8),%xmm14
   0x0000000000002fcc <+12236>:	mulps  0x3e0(%rcx),%xmm14
   0x0000000000002fe8 <+12264>:	addps  %xmm12,%xmm14
   0x0000000000003000 <+12288>:	movaps 0x7a0(%r8),%xmm14
   0x0000000000003008 <+12296>:	mulps  0x3e0(%rcx),%xmm14
   0x0000000000003038 <+12344>:	addps  %xmm11,%xmm14
   0x0000000000003050 <+12368>:	movaps 0x7e0(%r8),%xmm14
   0x0000000000003058 <+12376>:	mulps  0x3e0(%rcx),%xmm14
   0x0000000000003078 <+12408>:	addps  %xmm9,%xmm14

957
958	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
959	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
960	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002fb0 <+12208>:	movaps 0x730(%r8),%xmm10
   0x0000000000002fb8 <+12216>:	mulps  0x3f0(%rcx),%xmm10
   0x0000000000002fd4 <+12244>:	addps  %xmm13,%xmm10
   0x0000000000002fd8 <+12248>:	movaps 0x770(%r8),%xmm13
   0x0000000000002fe0 <+12256>:	mulps  0x3f0(%rcx),%xmm13
   0x0000000000002fec <+12268>:	movaps 0x7b0(%r8),%xmm12
   0x0000000000002ff4 <+12276>:	mulps  0x3f0(%rcx),%xmm12
   0x0000000000002ffc <+12284>:	addps  %xmm14,%xmm13
   0x000000000000303c <+12348>:	movaps 0x7f0(%r8),%xmm11
   0x0000000000003044 <+12356>:	mulps  0x3f0(%rcx),%xmm11
   0x000000000000304c <+12364>:	addps  %xmm14,%xmm12
   0x000000000000307c <+12412>:	addps  %xmm14,%xmm11

961	      }
962	    }
963
964	    /* Store C(2,4) block. */
965	    for (i = 0; i < 4; i++)
966	    {
967	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000003080 <+12416>:	mulps  %xmm0,%xmm10
   0x000000000000309c <+12444>:	mulps  %xmm0,%xmm13
   0x00000000000030ba <+12474>:	mulps  %xmm0,%xmm12
   0x00000000000030d2 <+12498>:	mulps  %xmm0,%xmm11

968	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_24]), C_row[i]);
   0x0000000000003084 <+12420>:	addps  0x1c0(%r9),%xmm10
   0x00000000000030a0 <+12448>:	addps  0x1d0(%r9),%xmm13
   0x00000000000030be <+12478>:	addps  0x1e0(%r9),%xmm12
   0x00000000000030d6 <+12502>:	addps  0x1f0(%r9),%xmm11

969	      _mm_store_ps(&C[i*4+C_OFFSET_24], C_row[i]);
   0x0000000000003090 <+12432>:	movaps %xmm10,0x1c0(%r9)
   0x00000000000030a8 <+12456>:	movaps %xmm13,0x1d0(%r9)
   0x00000000000030c6 <+12486>:	movaps %xmm12,0x1e0(%r9)
   0x00000000000030de <+12510>:	movaps %xmm11,0x1f0(%r9)

970	    }
971
972	    /* Reset C(3,1) matrix accumulators */
973	    C_row[0] = _mm_setzero_ps();
   0x00000000000030ce <+12494>:	xorps  %xmm12,%xmm12

974	    C_row[1] = _mm_setzero_ps();
   0x00000000000030e6 <+12518>:	xorps  %xmm11,%xmm11

975	    C_row[2] = _mm_setzero_ps();
   0x0000000000003098 <+12440>:	xorps  %xmm10,%xmm10

976	    C_row[3] = _mm_setzero_ps();
   0x000000000000308c <+12428>:	xorps  %xmm9,%xmm9

977
978	    if (norm_product[8][0] &&
   0x00000000000030b0 <+12464>:	movss  0x108(%rsp),%xmm13
   0x00000000000030ea <+12522>:	comiss %xmm1,%xmm13
   0x00000000000030ee <+12526>:	jb     0x35ec <stream_kernel+13804>

979	        norm_product[9][4] &&
   0x00000000000030f4 <+12532>:	movss  0xb0(%rsp),%xmm13
   0x00000000000030fe <+12542>:	comiss %xmm1,%xmm13
   0x0000000000003102 <+12546>:	jb     0x35ec <stream_kernel+13804>

980	        norm_product[10][8] &&
   0x0000000000003108 <+12552>:	movss  0xa8(%rsp),%xmm13
   0x0000000000003112 <+12562>:	comiss %xmm1,%xmm13
   0x0000000000003116 <+12566>:	jb     0x35ec <stream_kernel+13804>

981	        norm_product[11][12])
   0x000000000000311c <+12572>:	movss  0x98(%rsp),%xmm13
   0x0000000000003126 <+12582>:	comiss %xmm1,%xmm13
   0x000000000000312a <+12586>:	jb     0x35ec <stream_kernel+13804>

982	    {
983	      /* A(3,1)*B(1,1) = C(3,1). */
984	      for (i = 0; i < 4; i++)
985	      {
986	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
987	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
988	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003130 <+12592>:	movaps 0x800(%r8),%xmm14
   0x0000000000003150 <+12624>:	movaps 0x840(%r8),%xmm10
   0x0000000000003168 <+12648>:	mulps  (%rcx),%xmm14
   0x000000000000317b <+12667>:	mulps  (%rcx),%xmm10
   0x000000000000318d <+12685>:	movaps 0x8c0(%r8),%xmm14
   0x0000000000003195 <+12693>:	mulps  (%rcx),%xmm14
   0x00000000000031ae <+12718>:	movaps 0x880(%r8),%xmm11
   0x00000000000031b6 <+12726>:	mulps  (%rcx),%xmm11

989
990	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
991	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
992	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003138 <+12600>:	movaps 0x810(%r8),%xmm12
   0x0000000000003158 <+12632>:	movaps 0x850(%r8),%xmm9
   0x000000000000316c <+12652>:	mulps  0x10(%rcx),%xmm12
   0x000000000000317f <+12671>:	mulps  0x10(%rcx),%xmm9
   0x0000000000003189 <+12681>:	addps  %xmm14,%xmm12
   0x00000000000031ba <+12730>:	addps  %xmm10,%xmm9
   0x00000000000031be <+12734>:	movaps 0x890(%r8),%xmm10
   0x00000000000031c6 <+12742>:	mulps  0x10(%rcx),%xmm10
   0x00000000000031f0 <+12784>:	addps  %xmm11,%xmm10
   0x0000000000003205 <+12805>:	movaps 0x8d0(%r8),%xmm10
   0x000000000000320d <+12813>:	mulps  0x10(%rcx),%xmm10
   0x0000000000003223 <+12835>:	addps  %xmm14,%xmm10

993
994	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
995	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
996	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003140 <+12608>:	movaps 0x820(%r8),%xmm11
   0x0000000000003160 <+12640>:	movaps 0x860(%r8),%xmm15
   0x0000000000003171 <+12657>:	mulps  0x20(%rcx),%xmm11
   0x0000000000003184 <+12676>:	mulps  0x20(%rcx),%xmm15
   0x0000000000003199 <+12697>:	addps  %xmm12,%xmm11
   0x00000000000031cb <+12747>:	addps  %xmm9,%xmm15
   0x00000000000031cf <+12751>:	movaps 0x8a0(%r8),%xmm9
   0x00000000000031d7 <+12759>:	mulps  0x20(%rcx),%xmm9
   0x0000000000003201 <+12801>:	addps  %xmm10,%xmm9
   0x0000000000003216 <+12822>:	movaps 0x8e0(%r8),%xmm9
   0x000000000000321e <+12830>:	mulps  0x20(%rcx),%xmm9
   0x0000000000003237 <+12855>:	addps  %xmm10,%xmm9

997
998	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
999	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
1000	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003148 <+12616>:	movaps 0x830(%r8),%xmm13
   0x0000000000003176 <+12662>:	mulps  0x30(%rcx),%xmm13
   0x000000000000319d <+12701>:	movaps 0x870(%r8),%xmm12
   0x00000000000031a5 <+12709>:	mulps  0x30(%rcx),%xmm12
   0x00000000000031aa <+12714>:	addps  %xmm11,%xmm13
   0x00000000000031dc <+12764>:	addps  %xmm15,%xmm12
   0x00000000000031f4 <+12788>:	movaps 0x8b0(%r8),%xmm11
   0x00000000000031fc <+12796>:	mulps  0x30(%rcx),%xmm11
   0x0000000000003212 <+12818>:	addps  %xmm9,%xmm11
   0x000000000000323b <+12859>:	movaps 0x8f0(%r8),%xmm10
   0x0000000000003243 <+12867>:	mulps  0x30(%rcx),%xmm10
   0x000000000000325c <+12892>:	addps  %xmm9,%xmm10

1001	      }
1002
1003	      /* A(3,2)*B(2,1) = C(3,1). */
1004	      for (i = 0; i < 4; i++)
1005	      {
1006	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1007	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
1008	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000031e0 <+12768>:	movaps 0x900(%r8),%xmm15
   0x00000000000031e8 <+12776>:	mulps  0x100(%rcx),%xmm15
   0x0000000000003248 <+12872>:	addps  %xmm13,%xmm15
   0x0000000000003274 <+12916>:	movaps 0x940(%r8),%xmm15
   0x000000000000327c <+12924>:	mulps  0x100(%rcx),%xmm15
   0x000000000000329c <+12956>:	movaps 0x980(%r8),%xmm14
   0x00000000000032a4 <+12964>:	mulps  0x100(%rcx),%xmm14
   0x00000000000032ac <+12972>:	addps  %xmm12,%xmm15
   0x00000000000032c0 <+12992>:	addps  %xmm11,%xmm14
   0x0000000000003300 <+13056>:	movaps 0x9c0(%r8),%xmm13
   0x0000000000003308 <+13064>:	mulps  0x100(%rcx),%xmm13
   0x000000000000334c <+13132>:	addps  %xmm10,%xmm13

1009
1010	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1011	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
1012	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000324c <+12876>:	movaps 0x910(%r8),%xmm13
   0x0000000000003254 <+12884>:	mulps  0x110(%rcx),%xmm13
   0x0000000000003270 <+12912>:	addps  %xmm15,%xmm13
   0x00000000000032b0 <+12976>:	movaps 0x950(%r8),%xmm12
   0x00000000000032b8 <+12984>:	mulps  0x110(%rcx),%xmm12
   0x00000000000032c4 <+12996>:	movaps 0x990(%r8),%xmm11
   0x00000000000032cc <+13004>:	mulps  0x110(%rcx),%xmm11
   0x00000000000032d4 <+13012>:	addps  %xmm15,%xmm12
   0x0000000000003310 <+13072>:	addps  %xmm14,%xmm11
   0x0000000000003350 <+13136>:	movaps 0x9d0(%r8),%xmm10
   0x0000000000003358 <+13144>:	mulps  0x110(%rcx),%xmm10
   0x0000000000003374 <+13172>:	addps  %xmm13,%xmm10

1013
1014	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1015	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
1016	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003227 <+12839>:	movaps 0x920(%r8),%xmm14
   0x000000000000322f <+12847>:	mulps  0x120(%rcx),%xmm14
   0x0000000000003284 <+12932>:	addps  %xmm13,%xmm14
   0x0000000000003288 <+12936>:	movaps 0x960(%r8),%xmm13
   0x0000000000003290 <+12944>:	mulps  0x120(%rcx),%xmm13
   0x00000000000032e8 <+13032>:	addps  %xmm12,%xmm13
   0x00000000000032ec <+13036>:	movaps 0x9a0(%r8),%xmm12
   0x00000000000032f4 <+13044>:	mulps  0x120(%rcx),%xmm12
   0x0000000000003324 <+13092>:	addps  %xmm11,%xmm12
   0x000000000000333c <+13116>:	movaps 0x9e0(%r8),%xmm12
   0x0000000000003344 <+13124>:	mulps  0x120(%rcx),%xmm12
   0x0000000000003388 <+13192>:	addps  %xmm10,%xmm12

1017
1018	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1019	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
1020	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003260 <+12896>:	movaps 0x930(%r8),%xmm9
   0x0000000000003268 <+12904>:	mulps  0x130(%rcx),%xmm9
   0x0000000000003298 <+12952>:	addps  %xmm14,%xmm9
   0x00000000000032d8 <+13016>:	movaps 0x970(%r8),%xmm15
   0x00000000000032e0 <+13024>:	mulps  0x130(%rcx),%xmm15
   0x00000000000032fc <+13052>:	addps  %xmm13,%xmm15
   0x0000000000003314 <+13076>:	movaps 0x9b0(%r8),%xmm14
   0x000000000000331c <+13084>:	mulps  0x130(%rcx),%xmm14
   0x0000000000003338 <+13112>:	addps  %xmm12,%xmm14
   0x0000000000003378 <+13176>:	movaps 0x9f0(%r8),%xmm13
   0x0000000000003380 <+13184>:	mulps  0x130(%rcx),%xmm13
   0x000000000000339c <+13212>:	addps  %xmm12,%xmm13

1021	      }
1022
1023	      /* A(3,3)*B(3,1) = C(3,1). */
1024	      for (i = 0; i < 4; i++)
1025	      {
1026	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1027	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
1028	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003328 <+13096>:	movaps 0xa00(%r8),%xmm11
   0x0000000000003330 <+13104>:	mulps  0x200(%rcx),%xmm11
   0x0000000000003360 <+13152>:	addps  %xmm9,%xmm11
   0x000000000000338c <+13196>:	movaps 0xa40(%r8),%xmm10
   0x0000000000003394 <+13204>:	mulps  0x200(%rcx),%xmm10
   0x00000000000033c4 <+13252>:	addps  %xmm15,%xmm10
   0x00000000000033dc <+13276>:	movaps 0xa80(%r8),%xmm9
   0x00000000000033e4 <+13284>:	mulps  0x200(%rcx),%xmm9
   0x0000000000003414 <+13332>:	addps  %xmm14,%xmm9
   0x000000000000342c <+13356>:	movaps 0xac0(%r8),%xmm15
   0x0000000000003434 <+13364>:	mulps  0x200(%rcx),%xmm15
   0x0000000000003464 <+13412>:	addps  %xmm13,%xmm15

1029
1030	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1031	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
1032	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003364 <+13156>:	movaps 0xa10(%r8),%xmm9
   0x000000000000336c <+13164>:	mulps  0x210(%rcx),%xmm9
   0x00000000000033b0 <+13232>:	addps  %xmm11,%xmm9
   0x00000000000033c8 <+13256>:	movaps 0xa50(%r8),%xmm15
   0x00000000000033d0 <+13264>:	mulps  0x210(%rcx),%xmm15
   0x0000000000003400 <+13312>:	addps  %xmm10,%xmm15
   0x0000000000003418 <+13336>:	movaps 0xa90(%r8),%xmm14
   0x0000000000003420 <+13344>:	mulps  0x210(%rcx),%xmm14
   0x0000000000003450 <+13392>:	addps  %xmm9,%xmm14
   0x0000000000003468 <+13416>:	movaps 0xad0(%r8),%xmm13
   0x0000000000003470 <+13424>:	mulps  0x210(%rcx),%xmm13
   0x00000000000034a0 <+13472>:	addps  %xmm15,%xmm13

1033
1034	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1035	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
1036	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000033b4 <+13236>:	movaps 0xa20(%r8),%xmm11
   0x00000000000033bc <+13244>:	mulps  0x220(%rcx),%xmm11
   0x00000000000033d8 <+13272>:	addps  %xmm9,%xmm11
   0x0000000000003404 <+13316>:	movaps 0xa60(%r8),%xmm10
   0x000000000000340c <+13324>:	mulps  0x220(%rcx),%xmm10
   0x0000000000003428 <+13352>:	addps  %xmm15,%xmm10
   0x0000000000003454 <+13396>:	movaps 0xaa0(%r8),%xmm9
   0x000000000000345c <+13404>:	mulps  0x220(%rcx),%xmm9
   0x0000000000003478 <+13432>:	addps  %xmm14,%xmm9
   0x000000000000347c <+13436>:	movaps 0xae0(%r8),%xmm14
   0x0000000000003484 <+13444>:	mulps  0x220(%rcx),%xmm14
   0x00000000000034b4 <+13492>:	addps  %xmm13,%xmm14

1037
1038	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1039	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
1040	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000033a0 <+13216>:	movaps 0xa30(%r8),%xmm12
   0x00000000000033a8 <+13224>:	mulps  0x230(%rcx),%xmm12
   0x00000000000033ec <+13292>:	addps  %xmm11,%xmm12
   0x00000000000033f0 <+13296>:	movaps 0xa70(%r8),%xmm11
   0x00000000000033f8 <+13304>:	mulps  0x230(%rcx),%xmm11
   0x000000000000343c <+13372>:	addps  %xmm10,%xmm11
   0x0000000000003440 <+13376>:	movaps 0xab0(%r8),%xmm10
   0x0000000000003448 <+13384>:	mulps  0x230(%rcx),%xmm10
   0x000000000000348c <+13452>:	addps  %xmm9,%xmm10
   0x0000000000003490 <+13456>:	movaps 0xaf0(%r8),%xmm9
   0x0000000000003498 <+13464>:	mulps  0x230(%rcx),%xmm9
   0x00000000000034c8 <+13512>:	addps  %xmm14,%xmm9

1041	      }
1042
1043	      /* A(3,4)*B(4,1) = C(3,1). */
1044	      for (i = 0; i < 4; i++)
1045	      {
1046	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1047	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
1048	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000034a4 <+13476>:	movaps 0xb00(%r8),%xmm15
   0x00000000000034ac <+13484>:	mulps  0x300(%rcx),%xmm15
   0x00000000000034cc <+13516>:	movaps 0xb40(%r8),%xmm14
   0x00000000000034d4 <+13524>:	mulps  0x300(%rcx),%xmm14
   0x00000000000034dc <+13532>:	addps  %xmm12,%xmm15
   0x00000000000034f0 <+13552>:	addps  %xmm11,%xmm14
   0x0000000000003508 <+13576>:	movaps 0xb80(%r8),%xmm15
   0x0000000000003510 <+13584>:	mulps  0x300(%rcx),%xmm15
   0x0000000000003530 <+13616>:	movaps 0xbc0(%r8),%xmm14
   0x0000000000003538 <+13624>:	mulps  0x300(%rcx),%xmm14
   0x0000000000003554 <+13652>:	addps  %xmm10,%xmm15
   0x0000000000003580 <+13696>:	addps  %xmm9,%xmm14

1049
1050	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1051	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
1052	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000034e0 <+13536>:	movaps 0xb10(%r8),%xmm12
   0x00000000000034e8 <+13544>:	mulps  0x310(%rcx),%xmm12
   0x00000000000034f4 <+13556>:	movaps 0xb50(%r8),%xmm11
   0x00000000000034fc <+13564>:	mulps  0x310(%rcx),%xmm11
   0x0000000000003504 <+13572>:	addps  %xmm15,%xmm12
   0x000000000000352c <+13612>:	addps  %xmm14,%xmm11
   0x0000000000003558 <+13656>:	movaps 0xb90(%r8),%xmm10
   0x0000000000003560 <+13664>:	mulps  0x310(%rcx),%xmm10
   0x000000000000357c <+13692>:	addps  %xmm15,%xmm10
   0x0000000000003584 <+13700>:	movaps 0xbd0(%r8),%xmm9
   0x000000000000358c <+13708>:	mulps  0x310(%rcx),%xmm9
   0x00000000000035a8 <+13736>:	addps  %xmm14,%xmm9

1053
1054	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1055	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
1056	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000034b8 <+13496>:	movaps 0xb20(%r8),%xmm13
   0x00000000000034c0 <+13504>:	mulps  0x320(%rcx),%xmm13
   0x0000000000003518 <+13592>:	addps  %xmm12,%xmm13
   0x0000000000003544 <+13636>:	movaps 0xb60(%r8),%xmm13
   0x000000000000354c <+13644>:	mulps  0x320(%rcx),%xmm13
   0x0000000000003568 <+13672>:	addps  %xmm11,%xmm13
   0x0000000000003598 <+13720>:	movaps 0xba0(%r8),%xmm13
   0x00000000000035a0 <+13728>:	mulps  0x320(%rcx),%xmm13
   0x00000000000035ac <+13740>:	addps  %xmm10,%xmm13
   0x00000000000035c4 <+13764>:	movaps 0xbe0(%r8),%xmm13
   0x00000000000035cc <+13772>:	mulps  0x320(%rcx),%xmm13
   0x00000000000035d4 <+13780>:	addps  %xmm9,%xmm13

1057
1058	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1059	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
1060	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000351c <+13596>:	movaps 0xb30(%r8),%xmm12
   0x0000000000003524 <+13604>:	mulps  0x330(%rcx),%xmm12
   0x0000000000003540 <+13632>:	addps  %xmm13,%xmm12
   0x000000000000356c <+13676>:	movaps 0xb70(%r8),%xmm11
   0x0000000000003574 <+13684>:	mulps  0x330(%rcx),%xmm11
   0x0000000000003594 <+13716>:	addps  %xmm13,%xmm11
   0x00000000000035b0 <+13744>:	movaps 0xbb0(%r8),%xmm10
   0x00000000000035b8 <+13752>:	mulps  0x330(%rcx),%xmm10
   0x00000000000035c0 <+13760>:	addps  %xmm13,%xmm10
   0x00000000000035d8 <+13784>:	movaps 0xbf0(%r8),%xmm9
   0x00000000000035e0 <+13792>:	mulps  0x330(%rcx),%xmm9
   0x00000000000035e8 <+13800>:	addps  %xmm13,%xmm9

1061	      }
1062	    }
1063
1064	    /* Store C(3,1) block. */
1065	    for (i = 0; i < 4; i++)
1066	    {
1067	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000035ec <+13804>:	mulps  %xmm0,%xmm12
   0x000000000000360e <+13838>:	mulps  %xmm0,%xmm11
   0x0000000000003626 <+13862>:	mulps  %xmm0,%xmm10
   0x000000000000363e <+13886>:	mulps  %xmm0,%xmm9

1068	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_31]), C_row[i]);
   0x00000000000035f0 <+13808>:	addps  0x200(%r9),%xmm12
   0x0000000000003612 <+13842>:	addps  0x210(%r9),%xmm11
   0x000000000000362a <+13866>:	addps  0x220(%r9),%xmm10
   0x0000000000003642 <+13890>:	addps  0x230(%r9),%xmm9

1069	      _mm_store_ps(&C[i*4+C_OFFSET_31], C_row[i]);
   0x0000000000003602 <+13826>:	movaps %xmm12,0x200(%r9)
   0x000000000000361a <+13850>:	movaps %xmm11,0x210(%r9)
   0x0000000000003632 <+13874>:	movaps %xmm10,0x220(%r9)
   0x000000000000364a <+13898>:	movaps %xmm9,0x230(%r9)

1070	    }
1071
1072	    /* Reset C(3,2) matrix accumulators */
1073	    C_row[0] = _mm_setzero_ps();
   0x000000000000360a <+13834>:	xorps  %xmm12,%xmm12

1074	    C_row[1] = _mm_setzero_ps();
   0x0000000000003622 <+13858>:	xorps  %xmm11,%xmm11

1075	    C_row[2] = _mm_setzero_ps();
   0x000000000000363a <+13882>:	xorps  %xmm10,%xmm10

1076	    C_row[3] = _mm_setzero_ps();
   0x0000000000003652 <+13906>:	xorps  %xmm9,%xmm9

1077
1078	    if (norm_product[8][1] &&
   0x00000000000035f8 <+13816>:	movss  0xc0(%rsp),%xmm13
   0x0000000000003656 <+13910>:	comiss %xmm1,%xmm13
   0x000000000000365a <+13914>:	jb     0x3b59 <stream_kernel+15193>

1079	        norm_product[9][5] &&
   0x0000000000003660 <+13920>:	movss  0xa0(%rsp),%xmm13
   0x000000000000366a <+13930>:	comiss %xmm1,%xmm13
   0x000000000000366e <+13934>:	jb     0x3b59 <stream_kernel+15193>

1080	        norm_product[10][9] &&
   0x0000000000003674 <+13940>:	movss  0x88(%rsp),%xmm13
   0x000000000000367e <+13950>:	comiss %xmm1,%xmm13
   0x0000000000003682 <+13954>:	jb     0x3b59 <stream_kernel+15193>

1081	        norm_product[11][13])
   0x0000000000003688 <+13960>:	movss  0x78(%rsp),%xmm13
   0x000000000000368f <+13967>:	comiss %xmm1,%xmm13
   0x0000000000003693 <+13971>:	jb     0x3b59 <stream_kernel+15193>

1082	    {
1083	      /* A(3,1)*B(1,2) = C(3,2). */
1084	      for (i = 0; i < 4; i++)
1085	      {
1086	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
1087	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
1088	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003699 <+13977>:	movaps 0x800(%r8),%xmm14
   0x00000000000036b9 <+14009>:	movaps 0x840(%r8),%xmm10
   0x00000000000036d1 <+14033>:	mulps  0x40(%rcx),%xmm14
   0x00000000000036e5 <+14053>:	mulps  0x40(%rcx),%xmm10
   0x00000000000036f8 <+14072>:	movaps 0x8c0(%r8),%xmm14
   0x0000000000003700 <+14080>:	mulps  0x40(%rcx),%xmm14
   0x000000000000371a <+14106>:	movaps 0x880(%r8),%xmm11
   0x0000000000003722 <+14114>:	mulps  0x40(%rcx),%xmm11

1089
1090	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
1091	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
1092	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000036a1 <+13985>:	movaps 0x810(%r8),%xmm12
   0x00000000000036c1 <+14017>:	movaps 0x850(%r8),%xmm9
   0x00000000000036d6 <+14038>:	mulps  0x50(%rcx),%xmm12
   0x00000000000036ea <+14058>:	mulps  0x50(%rcx),%xmm9
   0x00000000000036f4 <+14068>:	addps  %xmm14,%xmm12
   0x0000000000003727 <+14119>:	addps  %xmm10,%xmm9
   0x000000000000372b <+14123>:	movaps 0x890(%r8),%xmm10
   0x0000000000003733 <+14131>:	mulps  0x50(%rcx),%xmm10
   0x000000000000375d <+14173>:	addps  %xmm11,%xmm10
   0x0000000000003772 <+14194>:	movaps 0x8d0(%r8),%xmm10
   0x000000000000377a <+14202>:	mulps  0x50(%rcx),%xmm10
   0x0000000000003790 <+14224>:	addps  %xmm14,%xmm10

1093
1094	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
1095	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
1096	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000036a9 <+13993>:	movaps 0x820(%r8),%xmm11
   0x00000000000036c9 <+14025>:	movaps 0x860(%r8),%xmm15
   0x00000000000036db <+14043>:	mulps  0x60(%rcx),%xmm11
   0x00000000000036ef <+14063>:	mulps  0x60(%rcx),%xmm15
   0x0000000000003705 <+14085>:	addps  %xmm12,%xmm11
   0x0000000000003738 <+14136>:	addps  %xmm9,%xmm15
   0x000000000000373c <+14140>:	movaps 0x8a0(%r8),%xmm9
   0x0000000000003744 <+14148>:	mulps  0x60(%rcx),%xmm9
   0x000000000000376e <+14190>:	addps  %xmm10,%xmm9
   0x0000000000003783 <+14211>:	movaps 0x8e0(%r8),%xmm9
   0x000000000000378b <+14219>:	mulps  0x60(%rcx),%xmm9
   0x00000000000037a4 <+14244>:	addps  %xmm10,%xmm9

1097
1098	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
1099	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
1100	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000036b1 <+14001>:	movaps 0x830(%r8),%xmm13
   0x00000000000036e0 <+14048>:	mulps  0x70(%rcx),%xmm13
   0x0000000000003709 <+14089>:	movaps 0x870(%r8),%xmm12
   0x0000000000003711 <+14097>:	mulps  0x70(%rcx),%xmm12
   0x0000000000003716 <+14102>:	addps  %xmm11,%xmm13
   0x0000000000003749 <+14153>:	addps  %xmm15,%xmm12
   0x0000000000003761 <+14177>:	movaps 0x8b0(%r8),%xmm11
   0x0000000000003769 <+14185>:	mulps  0x70(%rcx),%xmm11
   0x000000000000377f <+14207>:	addps  %xmm9,%xmm11
   0x00000000000037a8 <+14248>:	movaps 0x8f0(%r8),%xmm10
   0x00000000000037b0 <+14256>:	mulps  0x70(%rcx),%xmm10
   0x00000000000037c9 <+14281>:	addps  %xmm9,%xmm10

1101	      }
1102
1103	      /* A(3,2)*B(2,2) = C(3,2). */
1104	      for (i = 0; i < 4; i++)
1105	      {
1106	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1107	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
1108	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000374d <+14157>:	movaps 0x900(%r8),%xmm15
   0x0000000000003755 <+14165>:	mulps  0x140(%rcx),%xmm15
   0x00000000000037b5 <+14261>:	addps  %xmm13,%xmm15
   0x00000000000037e1 <+14305>:	movaps 0x940(%r8),%xmm15
   0x00000000000037e9 <+14313>:	mulps  0x140(%rcx),%xmm15
   0x0000000000003809 <+14345>:	movaps 0x980(%r8),%xmm14
   0x0000000000003811 <+14353>:	mulps  0x140(%rcx),%xmm14
   0x0000000000003819 <+14361>:	addps  %xmm12,%xmm15
   0x000000000000382d <+14381>:	addps  %xmm11,%xmm14
   0x000000000000386d <+14445>:	movaps 0x9c0(%r8),%xmm13
   0x0000000000003875 <+14453>:	mulps  0x140(%rcx),%xmm13
   0x00000000000038b9 <+14521>:	addps  %xmm10,%xmm13

1109
1110	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1111	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
1112	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000037b9 <+14265>:	movaps 0x910(%r8),%xmm13
   0x00000000000037c1 <+14273>:	mulps  0x150(%rcx),%xmm13
   0x00000000000037dd <+14301>:	addps  %xmm15,%xmm13
   0x000000000000381d <+14365>:	movaps 0x950(%r8),%xmm12
   0x0000000000003825 <+14373>:	mulps  0x150(%rcx),%xmm12
   0x0000000000003831 <+14385>:	movaps 0x990(%r8),%xmm11
   0x0000000000003839 <+14393>:	mulps  0x150(%rcx),%xmm11
   0x0000000000003841 <+14401>:	addps  %xmm15,%xmm12
   0x000000000000387d <+14461>:	addps  %xmm14,%xmm11
   0x00000000000038bd <+14525>:	movaps 0x9d0(%r8),%xmm10
   0x00000000000038c5 <+14533>:	mulps  0x150(%rcx),%xmm10
   0x00000000000038e1 <+14561>:	addps  %xmm13,%xmm10

1113
1114	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1115	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
1116	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003794 <+14228>:	movaps 0x920(%r8),%xmm14
   0x000000000000379c <+14236>:	mulps  0x160(%rcx),%xmm14
   0x00000000000037f1 <+14321>:	addps  %xmm13,%xmm14
   0x00000000000037f5 <+14325>:	movaps 0x960(%r8),%xmm13
   0x00000000000037fd <+14333>:	mulps  0x160(%rcx),%xmm13
   0x0000000000003855 <+14421>:	addps  %xmm12,%xmm13
   0x0000000000003859 <+14425>:	movaps 0x9a0(%r8),%xmm12
   0x0000000000003861 <+14433>:	mulps  0x160(%rcx),%xmm12
   0x0000000000003891 <+14481>:	addps  %xmm11,%xmm12
   0x00000000000038a9 <+14505>:	movaps 0x9e0(%r8),%xmm12
   0x00000000000038b1 <+14513>:	mulps  0x160(%rcx),%xmm12
   0x00000000000038f5 <+14581>:	addps  %xmm10,%xmm12

1117
1118	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1119	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
1120	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000037cd <+14285>:	movaps 0x930(%r8),%xmm9
   0x00000000000037d5 <+14293>:	mulps  0x170(%rcx),%xmm9
   0x0000000000003805 <+14341>:	addps  %xmm14,%xmm9
   0x0000000000003845 <+14405>:	movaps 0x970(%r8),%xmm15
   0x000000000000384d <+14413>:	mulps  0x170(%rcx),%xmm15
   0x0000000000003869 <+14441>:	addps  %xmm13,%xmm15
   0x0000000000003881 <+14465>:	movaps 0x9b0(%r8),%xmm14
   0x0000000000003889 <+14473>:	mulps  0x170(%rcx),%xmm14
   0x00000000000038a5 <+14501>:	addps  %xmm12,%xmm14
   0x00000000000038e5 <+14565>:	movaps 0x9f0(%r8),%xmm13
   0x00000000000038ed <+14573>:	mulps  0x170(%rcx),%xmm13
   0x0000000000003909 <+14601>:	addps  %xmm12,%xmm13

1121	      }
1122
1123	      /* A(3,3)*B(3,2) = C(3,2). */
1124	      for (i = 0; i < 4; i++)
1125	      {
1126	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1127	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
1128	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003895 <+14485>:	movaps 0xa00(%r8),%xmm11
   0x000000000000389d <+14493>:	mulps  0x240(%rcx),%xmm11
   0x00000000000038cd <+14541>:	addps  %xmm9,%xmm11
   0x00000000000038f9 <+14585>:	movaps 0xa40(%r8),%xmm10
   0x0000000000003901 <+14593>:	mulps  0x240(%rcx),%xmm10
   0x0000000000003931 <+14641>:	addps  %xmm15,%xmm10
   0x0000000000003949 <+14665>:	movaps 0xa80(%r8),%xmm9
   0x0000000000003951 <+14673>:	mulps  0x240(%rcx),%xmm9
   0x0000000000003981 <+14721>:	addps  %xmm14,%xmm9
   0x0000000000003999 <+14745>:	movaps 0xac0(%r8),%xmm15
   0x00000000000039a1 <+14753>:	mulps  0x240(%rcx),%xmm15
   0x00000000000039d1 <+14801>:	addps  %xmm13,%xmm15

1129
1130	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1131	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
1132	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000038d1 <+14545>:	movaps 0xa10(%r8),%xmm9
   0x00000000000038d9 <+14553>:	mulps  0x250(%rcx),%xmm9
   0x000000000000391d <+14621>:	addps  %xmm11,%xmm9
   0x0000000000003935 <+14645>:	movaps 0xa50(%r8),%xmm15
   0x000000000000393d <+14653>:	mulps  0x250(%rcx),%xmm15
   0x000000000000396d <+14701>:	addps  %xmm10,%xmm15
   0x0000000000003985 <+14725>:	movaps 0xa90(%r8),%xmm14
   0x000000000000398d <+14733>:	mulps  0x250(%rcx),%xmm14
   0x00000000000039bd <+14781>:	addps  %xmm9,%xmm14
   0x00000000000039d5 <+14805>:	movaps 0xad0(%r8),%xmm13
   0x00000000000039dd <+14813>:	mulps  0x250(%rcx),%xmm13
   0x0000000000003a0d <+14861>:	addps  %xmm15,%xmm13

1133
1134	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1135	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
1136	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003921 <+14625>:	movaps 0xa20(%r8),%xmm11
   0x0000000000003929 <+14633>:	mulps  0x260(%rcx),%xmm11
   0x0000000000003945 <+14661>:	addps  %xmm9,%xmm11
   0x0000000000003971 <+14705>:	movaps 0xa60(%r8),%xmm10
   0x0000000000003979 <+14713>:	mulps  0x260(%rcx),%xmm10
   0x0000000000003995 <+14741>:	addps  %xmm15,%xmm10
   0x00000000000039c1 <+14785>:	movaps 0xaa0(%r8),%xmm9
   0x00000000000039c9 <+14793>:	mulps  0x260(%rcx),%xmm9
   0x00000000000039e5 <+14821>:	addps  %xmm14,%xmm9
   0x00000000000039e9 <+14825>:	movaps 0xae0(%r8),%xmm14
   0x00000000000039f1 <+14833>:	mulps  0x260(%rcx),%xmm14
   0x0000000000003a21 <+14881>:	addps  %xmm13,%xmm14

1137
1138	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1139	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
1140	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000390d <+14605>:	movaps 0xa30(%r8),%xmm12
   0x0000000000003915 <+14613>:	mulps  0x270(%rcx),%xmm12
   0x0000000000003959 <+14681>:	addps  %xmm11,%xmm12
   0x000000000000395d <+14685>:	movaps 0xa70(%r8),%xmm11
   0x0000000000003965 <+14693>:	mulps  0x270(%rcx),%xmm11
   0x00000000000039a9 <+14761>:	addps  %xmm10,%xmm11
   0x00000000000039ad <+14765>:	movaps 0xab0(%r8),%xmm10
   0x00000000000039b5 <+14773>:	mulps  0x270(%rcx),%xmm10
   0x00000000000039f9 <+14841>:	addps  %xmm9,%xmm10
   0x00000000000039fd <+14845>:	movaps 0xaf0(%r8),%xmm9
   0x0000000000003a05 <+14853>:	mulps  0x270(%rcx),%xmm9
   0x0000000000003a35 <+14901>:	addps  %xmm14,%xmm9

1141	      }
1142
1143	      /* A(3,4)*B(4,2) = C(3,2). */
1144	      for (i = 0; i < 4; i++)
1145	      {
1146	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1147	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
1148	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003a11 <+14865>:	movaps 0xb00(%r8),%xmm15
   0x0000000000003a19 <+14873>:	mulps  0x340(%rcx),%xmm15
   0x0000000000003a39 <+14905>:	movaps 0xb40(%r8),%xmm14
   0x0000000000003a41 <+14913>:	mulps  0x340(%rcx),%xmm14
   0x0000000000003a49 <+14921>:	addps  %xmm12,%xmm15
   0x0000000000003a5d <+14941>:	addps  %xmm11,%xmm14
   0x0000000000003a75 <+14965>:	movaps 0xb80(%r8),%xmm15
   0x0000000000003a7d <+14973>:	mulps  0x340(%rcx),%xmm15
   0x0000000000003a9d <+15005>:	movaps 0xbc0(%r8),%xmm14
   0x0000000000003aa5 <+15013>:	mulps  0x340(%rcx),%xmm14
   0x0000000000003ac1 <+15041>:	addps  %xmm10,%xmm15
   0x0000000000003aed <+15085>:	addps  %xmm9,%xmm14

1149
1150	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1151	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
1152	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003a4d <+14925>:	movaps 0xb10(%r8),%xmm12
   0x0000000000003a55 <+14933>:	mulps  0x350(%rcx),%xmm12
   0x0000000000003a61 <+14945>:	movaps 0xb50(%r8),%xmm11
   0x0000000000003a69 <+14953>:	mulps  0x350(%rcx),%xmm11
   0x0000000000003a71 <+14961>:	addps  %xmm15,%xmm12
   0x0000000000003a99 <+15001>:	addps  %xmm14,%xmm11
   0x0000000000003ac5 <+15045>:	movaps 0xb90(%r8),%xmm10
   0x0000000000003acd <+15053>:	mulps  0x350(%rcx),%xmm10
   0x0000000000003ae9 <+15081>:	addps  %xmm15,%xmm10
   0x0000000000003af1 <+15089>:	movaps 0xbd0(%r8),%xmm9
   0x0000000000003af9 <+15097>:	mulps  0x350(%rcx),%xmm9
   0x0000000000003b15 <+15125>:	addps  %xmm14,%xmm9

1153
1154	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1155	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
1156	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003a25 <+14885>:	movaps 0xb20(%r8),%xmm13
   0x0000000000003a2d <+14893>:	mulps  0x360(%rcx),%xmm13
   0x0000000000003a85 <+14981>:	addps  %xmm12,%xmm13
   0x0000000000003ab1 <+15025>:	movaps 0xb60(%r8),%xmm13
   0x0000000000003ab9 <+15033>:	mulps  0x360(%rcx),%xmm13
   0x0000000000003ad5 <+15061>:	addps  %xmm11,%xmm13
   0x0000000000003b05 <+15109>:	movaps 0xba0(%r8),%xmm13
   0x0000000000003b0d <+15117>:	mulps  0x360(%rcx),%xmm13
   0x0000000000003b19 <+15129>:	addps  %xmm10,%xmm13
   0x0000000000003b31 <+15153>:	movaps 0xbe0(%r8),%xmm13
   0x0000000000003b39 <+15161>:	mulps  0x360(%rcx),%xmm13
   0x0000000000003b41 <+15169>:	addps  %xmm9,%xmm13

1157
1158	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1159	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
1160	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003a89 <+14985>:	movaps 0xb30(%r8),%xmm12
   0x0000000000003a91 <+14993>:	mulps  0x370(%rcx),%xmm12
   0x0000000000003aad <+15021>:	addps  %xmm13,%xmm12
   0x0000000000003ad9 <+15065>:	movaps 0xb70(%r8),%xmm11
   0x0000000000003ae1 <+15073>:	mulps  0x370(%rcx),%xmm11
   0x0000000000003b01 <+15105>:	addps  %xmm13,%xmm11
   0x0000000000003b1d <+15133>:	movaps 0xbb0(%r8),%xmm10
   0x0000000000003b25 <+15141>:	mulps  0x370(%rcx),%xmm10
   0x0000000000003b2d <+15149>:	addps  %xmm13,%xmm10
   0x0000000000003b45 <+15173>:	movaps 0xbf0(%r8),%xmm9
   0x0000000000003b4d <+15181>:	mulps  0x370(%rcx),%xmm9
   0x0000000000003b55 <+15189>:	addps  %xmm13,%xmm9

1161	      }
1162	    }
1163
1164	    /* Store C(3,2) block. */
1165	    for (i = 0; i < 4; i++)
1166	    {
1167	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000003b59 <+15193>:	mulps  %xmm0,%xmm12
   0x0000000000003b7b <+15227>:	mulps  %xmm0,%xmm11
   0x0000000000003b93 <+15251>:	mulps  %xmm0,%xmm10
   0x0000000000003bab <+15275>:	mulps  %xmm0,%xmm9

1168	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_32]), C_row[i]);
   0x0000000000003b5d <+15197>:	addps  0x240(%r9),%xmm12
   0x0000000000003b7f <+15231>:	addps  0x250(%r9),%xmm11
   0x0000000000003b97 <+15255>:	addps  0x260(%r9),%xmm10
   0x0000000000003baf <+15279>:	addps  0x270(%r9),%xmm9

1169	      _mm_store_ps(&C[i*4+C_OFFSET_32], C_row[i]);
   0x0000000000003b6f <+15215>:	movaps %xmm12,0x240(%r9)
   0x0000000000003b87 <+15239>:	movaps %xmm11,0x250(%r9)
   0x0000000000003b9f <+15263>:	movaps %xmm10,0x260(%r9)
   0x0000000000003bb7 <+15287>:	movaps %xmm9,0x270(%r9)

1170	    }
1171
1172	    /* Reset C(3,3) matrix accumulators */
1173	    C_row[0] = _mm_setzero_ps();
   0x0000000000003b77 <+15223>:	xorps  %xmm12,%xmm12

1174	    C_row[1] = _mm_setzero_ps();
   0x0000000000003b8f <+15247>:	xorps  %xmm11,%xmm11

1175	    C_row[2] = _mm_setzero_ps();
   0x0000000000003ba7 <+15271>:	xorps  %xmm10,%xmm10

1176	    C_row[3] = _mm_setzero_ps();
   0x0000000000003bbf <+15295>:	xorps  %xmm9,%xmm9

1177
1178	    if (norm_product[8][2] &&
   0x0000000000003b65 <+15205>:	movss  0x128(%rsp),%xmm13
   0x0000000000003bc3 <+15299>:	comiss %xmm1,%xmm13
   0x0000000000003bc7 <+15303>:	jb     0x40f3 <stream_kernel+16627>

1179	        norm_product[9][6] &&
   0x0000000000003bcd <+15309>:	movss  0x80(%rsp),%xmm13
   0x0000000000003bd7 <+15319>:	comiss %xmm1,%xmm13
   0x0000000000003bdb <+15323>:	jb     0x40f3 <stream_kernel+16627>

1180	        norm_product[10][10] &&
   0x0000000000003be1 <+15329>:	movss  0x60(%rsp),%xmm13
   0x0000000000003be8 <+15336>:	comiss %xmm1,%xmm13
   0x0000000000003bec <+15340>:	jb     0x40f3 <stream_kernel+16627>

1181	        norm_product[11][14])
   0x0000000000003bf2 <+15346>:	movss  0x50(%rsp),%xmm13
   0x0000000000003bf9 <+15353>:	comiss %xmm1,%xmm13
   0x0000000000003bfd <+15357>:	jb     0x40f3 <stream_kernel+16627>

1182	    {
1183	      /* A(3,1)*B(1,3) = C(3,3). */
1184	      for (i = 0; i < 4; i++)
1185	      {
1186	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
1187	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
1188	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003c03 <+15363>:	movaps 0x800(%r8),%xmm14
   0x0000000000003c23 <+15395>:	movaps 0x840(%r8),%xmm10
   0x0000000000003c3b <+15419>:	mulps  0x80(%rcx),%xmm14
   0x0000000000003c5b <+15451>:	mulps  0x80(%rcx),%xmm10
   0x0000000000003c77 <+15479>:	movaps 0x8c0(%r8),%xmm14
   0x0000000000003c7f <+15487>:	mulps  0x80(%rcx),%xmm14
   0x0000000000003c9f <+15519>:	movaps 0x880(%r8),%xmm11
   0x0000000000003ca7 <+15527>:	mulps  0x80(%rcx),%xmm11

1189
1190	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
1191	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
1192	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003c0b <+15371>:	movaps 0x810(%r8),%xmm12
   0x0000000000003c2b <+15403>:	movaps 0x850(%r8),%xmm9
   0x0000000000003c43 <+15427>:	mulps  0x90(%rcx),%xmm12
   0x0000000000003c63 <+15459>:	mulps  0x90(%rcx),%xmm9
   0x0000000000003c73 <+15475>:	addps  %xmm14,%xmm12
   0x0000000000003caf <+15535>:	addps  %xmm10,%xmm9
   0x0000000000003cb3 <+15539>:	movaps 0x890(%r8),%xmm10
   0x0000000000003cbb <+15547>:	mulps  0x90(%rcx),%xmm10
   0x0000000000003ceb <+15595>:	addps  %xmm11,%xmm10
   0x0000000000003d03 <+15619>:	movaps 0x8d0(%r8),%xmm10
   0x0000000000003d0b <+15627>:	mulps  0x90(%rcx),%xmm10
   0x0000000000003d27 <+15655>:	addps  %xmm14,%xmm10

1193
1194	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
1195	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
1196	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003c13 <+15379>:	movaps 0x820(%r8),%xmm11
   0x0000000000003c33 <+15411>:	movaps 0x860(%r8),%xmm15
   0x0000000000003c4b <+15435>:	mulps  0xa0(%rcx),%xmm11
   0x0000000000003c6b <+15467>:	mulps  0xa0(%rcx),%xmm15
   0x0000000000003c87 <+15495>:	addps  %xmm12,%xmm11
   0x0000000000003cc3 <+15555>:	addps  %xmm9,%xmm15
   0x0000000000003cc7 <+15559>:	movaps 0x8a0(%r8),%xmm9
   0x0000000000003ccf <+15567>:	mulps  0xa0(%rcx),%xmm9
   0x0000000000003cff <+15615>:	addps  %xmm10,%xmm9
   0x0000000000003d17 <+15639>:	movaps 0x8e0(%r8),%xmm9
   0x0000000000003d1f <+15647>:	mulps  0xa0(%rcx),%xmm9
   0x0000000000003d3b <+15675>:	addps  %xmm10,%xmm9

1197
1198	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
1199	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
1200	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003c1b <+15387>:	movaps 0x830(%r8),%xmm13
   0x0000000000003c53 <+15443>:	mulps  0xb0(%rcx),%xmm13
   0x0000000000003c8b <+15499>:	movaps 0x870(%r8),%xmm12
   0x0000000000003c93 <+15507>:	mulps  0xb0(%rcx),%xmm12
   0x0000000000003c9b <+15515>:	addps  %xmm11,%xmm13
   0x0000000000003cd7 <+15575>:	addps  %xmm15,%xmm12
   0x0000000000003cef <+15599>:	movaps 0x8b0(%r8),%xmm11
   0x0000000000003cf7 <+15607>:	mulps  0xb0(%rcx),%xmm11
   0x0000000000003d13 <+15635>:	addps  %xmm9,%xmm11
   0x0000000000003d3f <+15679>:	movaps 0x8f0(%r8),%xmm10
   0x0000000000003d47 <+15687>:	mulps  0xb0(%rcx),%xmm10
   0x0000000000003d63 <+15715>:	addps  %xmm9,%xmm10

1201	      }
1202
1203	      /* A(3,2)*B(2,3) = C(3,3). */
1204	      for (i = 0; i < 4; i++)
1205	      {
1206	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1207	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
1208	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003cdb <+15579>:	movaps 0x900(%r8),%xmm15
   0x0000000000003ce3 <+15587>:	mulps  0x180(%rcx),%xmm15
   0x0000000000003d4f <+15695>:	addps  %xmm13,%xmm15
   0x0000000000003d7b <+15739>:	movaps 0x940(%r8),%xmm15
   0x0000000000003d83 <+15747>:	mulps  0x180(%rcx),%xmm15
   0x0000000000003da3 <+15779>:	movaps 0x980(%r8),%xmm14
   0x0000000000003dab <+15787>:	mulps  0x180(%rcx),%xmm14
   0x0000000000003db3 <+15795>:	addps  %xmm12,%xmm15
   0x0000000000003dc7 <+15815>:	addps  %xmm11,%xmm14
   0x0000000000003e07 <+15879>:	movaps 0x9c0(%r8),%xmm13
   0x0000000000003e0f <+15887>:	mulps  0x180(%rcx),%xmm13
   0x0000000000003e53 <+15955>:	addps  %xmm10,%xmm13

1209
1210	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1211	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
1212	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003d53 <+15699>:	movaps 0x910(%r8),%xmm13
   0x0000000000003d5b <+15707>:	mulps  0x190(%rcx),%xmm13
   0x0000000000003d77 <+15735>:	addps  %xmm15,%xmm13
   0x0000000000003db7 <+15799>:	movaps 0x950(%r8),%xmm12
   0x0000000000003dbf <+15807>:	mulps  0x190(%rcx),%xmm12
   0x0000000000003dcb <+15819>:	movaps 0x990(%r8),%xmm11
   0x0000000000003dd3 <+15827>:	mulps  0x190(%rcx),%xmm11
   0x0000000000003ddb <+15835>:	addps  %xmm15,%xmm12
   0x0000000000003e17 <+15895>:	addps  %xmm14,%xmm11
   0x0000000000003e57 <+15959>:	movaps 0x9d0(%r8),%xmm10
   0x0000000000003e5f <+15967>:	mulps  0x190(%rcx),%xmm10
   0x0000000000003e7b <+15995>:	addps  %xmm13,%xmm10

1213
1214	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1215	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
1216	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003d2b <+15659>:	movaps 0x920(%r8),%xmm14
   0x0000000000003d33 <+15667>:	mulps  0x1a0(%rcx),%xmm14
   0x0000000000003d8b <+15755>:	addps  %xmm13,%xmm14
   0x0000000000003d8f <+15759>:	movaps 0x960(%r8),%xmm13
   0x0000000000003d97 <+15767>:	mulps  0x1a0(%rcx),%xmm13
   0x0000000000003def <+15855>:	addps  %xmm12,%xmm13
   0x0000000000003df3 <+15859>:	movaps 0x9a0(%r8),%xmm12
   0x0000000000003dfb <+15867>:	mulps  0x1a0(%rcx),%xmm12
   0x0000000000003e2b <+15915>:	addps  %xmm11,%xmm12
   0x0000000000003e43 <+15939>:	movaps 0x9e0(%r8),%xmm12
   0x0000000000003e4b <+15947>:	mulps  0x1a0(%rcx),%xmm12
   0x0000000000003e8f <+16015>:	addps  %xmm10,%xmm12

1217
1218	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1219	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
1220	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003d67 <+15719>:	movaps 0x930(%r8),%xmm9
   0x0000000000003d6f <+15727>:	mulps  0x1b0(%rcx),%xmm9
   0x0000000000003d9f <+15775>:	addps  %xmm14,%xmm9
   0x0000000000003ddf <+15839>:	movaps 0x970(%r8),%xmm15
   0x0000000000003de7 <+15847>:	mulps  0x1b0(%rcx),%xmm15
   0x0000000000003e03 <+15875>:	addps  %xmm13,%xmm15
   0x0000000000003e1b <+15899>:	movaps 0x9b0(%r8),%xmm14
   0x0000000000003e23 <+15907>:	mulps  0x1b0(%rcx),%xmm14
   0x0000000000003e3f <+15935>:	addps  %xmm12,%xmm14
   0x0000000000003e7f <+15999>:	movaps 0x9f0(%r8),%xmm13
   0x0000000000003e87 <+16007>:	mulps  0x1b0(%rcx),%xmm13
   0x0000000000003ea3 <+16035>:	addps  %xmm12,%xmm13

1221	      }
1222
1223	      /* A(3,3)*B(3,3) = C(3,3). */
1224	      for (i = 0; i < 4; i++)
1225	      {
1226	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1227	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
1228	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003e2f <+15919>:	movaps 0xa00(%r8),%xmm11
   0x0000000000003e37 <+15927>:	mulps  0x280(%rcx),%xmm11
   0x0000000000003e67 <+15975>:	addps  %xmm9,%xmm11
   0x0000000000003e93 <+16019>:	movaps 0xa40(%r8),%xmm10
   0x0000000000003e9b <+16027>:	mulps  0x280(%rcx),%xmm10
   0x0000000000003ecb <+16075>:	addps  %xmm15,%xmm10
   0x0000000000003ee3 <+16099>:	movaps 0xa80(%r8),%xmm9
   0x0000000000003eeb <+16107>:	mulps  0x280(%rcx),%xmm9
   0x0000000000003f1b <+16155>:	addps  %xmm14,%xmm9
   0x0000000000003f33 <+16179>:	movaps 0xac0(%r8),%xmm15
   0x0000000000003f3b <+16187>:	mulps  0x280(%rcx),%xmm15
   0x0000000000003f6b <+16235>:	addps  %xmm13,%xmm15

1229
1230	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1231	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
1232	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003e6b <+15979>:	movaps 0xa10(%r8),%xmm9
   0x0000000000003e73 <+15987>:	mulps  0x290(%rcx),%xmm9
   0x0000000000003eb7 <+16055>:	addps  %xmm11,%xmm9
   0x0000000000003ecf <+16079>:	movaps 0xa50(%r8),%xmm15
   0x0000000000003ed7 <+16087>:	mulps  0x290(%rcx),%xmm15
   0x0000000000003f07 <+16135>:	addps  %xmm10,%xmm15
   0x0000000000003f1f <+16159>:	movaps 0xa90(%r8),%xmm14
   0x0000000000003f27 <+16167>:	mulps  0x290(%rcx),%xmm14
   0x0000000000003f57 <+16215>:	addps  %xmm9,%xmm14
   0x0000000000003f6f <+16239>:	movaps 0xad0(%r8),%xmm13
   0x0000000000003f77 <+16247>:	mulps  0x290(%rcx),%xmm13
   0x0000000000003fa7 <+16295>:	addps  %xmm15,%xmm13

1233
1234	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1235	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
1236	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003ebb <+16059>:	movaps 0xa20(%r8),%xmm11
   0x0000000000003ec3 <+16067>:	mulps  0x2a0(%rcx),%xmm11
   0x0000000000003edf <+16095>:	addps  %xmm9,%xmm11
   0x0000000000003f0b <+16139>:	movaps 0xa60(%r8),%xmm10
   0x0000000000003f13 <+16147>:	mulps  0x2a0(%rcx),%xmm10
   0x0000000000003f2f <+16175>:	addps  %xmm15,%xmm10
   0x0000000000003f5b <+16219>:	movaps 0xaa0(%r8),%xmm9
   0x0000000000003f63 <+16227>:	mulps  0x2a0(%rcx),%xmm9
   0x0000000000003f7f <+16255>:	addps  %xmm14,%xmm9
   0x0000000000003f83 <+16259>:	movaps 0xae0(%r8),%xmm14
   0x0000000000003f8b <+16267>:	mulps  0x2a0(%rcx),%xmm14
   0x0000000000003fbb <+16315>:	addps  %xmm13,%xmm14

1237
1238	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1239	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
1240	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003ea7 <+16039>:	movaps 0xa30(%r8),%xmm12
   0x0000000000003eaf <+16047>:	mulps  0x2b0(%rcx),%xmm12
   0x0000000000003ef3 <+16115>:	addps  %xmm11,%xmm12
   0x0000000000003ef7 <+16119>:	movaps 0xa70(%r8),%xmm11
   0x0000000000003eff <+16127>:	mulps  0x2b0(%rcx),%xmm11
   0x0000000000003f43 <+16195>:	addps  %xmm10,%xmm11
   0x0000000000003f47 <+16199>:	movaps 0xab0(%r8),%xmm10
   0x0000000000003f4f <+16207>:	mulps  0x2b0(%rcx),%xmm10
   0x0000000000003f93 <+16275>:	addps  %xmm9,%xmm10
   0x0000000000003f97 <+16279>:	movaps 0xaf0(%r8),%xmm9
   0x0000000000003f9f <+16287>:	mulps  0x2b0(%rcx),%xmm9
   0x0000000000003fcf <+16335>:	addps  %xmm14,%xmm9

1241	      }
1242
1243	      /* A(3,4)*B(4,3) = C(3,3). */
1244	      for (i = 0; i < 4; i++)
1245	      {
1246	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1247	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
1248	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003fab <+16299>:	movaps 0xb00(%r8),%xmm15
   0x0000000000003fb3 <+16307>:	mulps  0x380(%rcx),%xmm15
   0x0000000000003fd3 <+16339>:	movaps 0xb40(%r8),%xmm14
   0x0000000000003fdb <+16347>:	mulps  0x380(%rcx),%xmm14
   0x0000000000003fe3 <+16355>:	addps  %xmm12,%xmm15
   0x0000000000003ff7 <+16375>:	addps  %xmm11,%xmm14
   0x000000000000400f <+16399>:	movaps 0xb80(%r8),%xmm15
   0x0000000000004017 <+16407>:	mulps  0x380(%rcx),%xmm15
   0x0000000000004037 <+16439>:	movaps 0xbc0(%r8),%xmm14
   0x000000000000403f <+16447>:	mulps  0x380(%rcx),%xmm14
   0x000000000000405b <+16475>:	addps  %xmm10,%xmm15
   0x0000000000004087 <+16519>:	addps  %xmm9,%xmm14

1249
1250	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1251	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
1252	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003fe7 <+16359>:	movaps 0xb10(%r8),%xmm12
   0x0000000000003fef <+16367>:	mulps  0x390(%rcx),%xmm12
   0x0000000000003ffb <+16379>:	movaps 0xb50(%r8),%xmm11
   0x0000000000004003 <+16387>:	mulps  0x390(%rcx),%xmm11
   0x000000000000400b <+16395>:	addps  %xmm15,%xmm12
   0x0000000000004033 <+16435>:	addps  %xmm14,%xmm11
   0x000000000000405f <+16479>:	movaps 0xb90(%r8),%xmm10
   0x0000000000004067 <+16487>:	mulps  0x390(%rcx),%xmm10
   0x0000000000004083 <+16515>:	addps  %xmm15,%xmm10
   0x000000000000408b <+16523>:	movaps 0xbd0(%r8),%xmm9
   0x0000000000004093 <+16531>:	mulps  0x390(%rcx),%xmm9
   0x00000000000040af <+16559>:	addps  %xmm14,%xmm9

1253
1254	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1255	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
1256	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003fbf <+16319>:	movaps 0xb20(%r8),%xmm13
   0x0000000000003fc7 <+16327>:	mulps  0x3a0(%rcx),%xmm13
   0x000000000000401f <+16415>:	addps  %xmm12,%xmm13
   0x000000000000404b <+16459>:	movaps 0xb60(%r8),%xmm13
   0x0000000000004053 <+16467>:	mulps  0x3a0(%rcx),%xmm13
   0x000000000000406f <+16495>:	addps  %xmm11,%xmm13
   0x000000000000409f <+16543>:	movaps 0xba0(%r8),%xmm13
   0x00000000000040a7 <+16551>:	mulps  0x3a0(%rcx),%xmm13
   0x00000000000040b3 <+16563>:	addps  %xmm10,%xmm13
   0x00000000000040cb <+16587>:	movaps 0xbe0(%r8),%xmm13
   0x00000000000040d3 <+16595>:	mulps  0x3a0(%rcx),%xmm13
   0x00000000000040db <+16603>:	addps  %xmm9,%xmm13

1257
1258	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1259	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
1260	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004023 <+16419>:	movaps 0xb30(%r8),%xmm12
   0x000000000000402b <+16427>:	mulps  0x3b0(%rcx),%xmm12
   0x0000000000004047 <+16455>:	addps  %xmm13,%xmm12
   0x0000000000004073 <+16499>:	movaps 0xb70(%r8),%xmm11
   0x000000000000407b <+16507>:	mulps  0x3b0(%rcx),%xmm11
   0x000000000000409b <+16539>:	addps  %xmm13,%xmm11
   0x00000000000040b7 <+16567>:	movaps 0xbb0(%r8),%xmm10
   0x00000000000040bf <+16575>:	mulps  0x3b0(%rcx),%xmm10
   0x00000000000040c7 <+16583>:	addps  %xmm13,%xmm10
   0x00000000000040df <+16607>:	movaps 0xbf0(%r8),%xmm9
   0x00000000000040e7 <+16615>:	mulps  0x3b0(%rcx),%xmm9
   0x00000000000040ef <+16623>:	addps  %xmm13,%xmm9

1261	      }
1262	    }
1263
1264	    /* Store C(3,3) block. */
1265	    for (i = 0; i < 4; i++)
1266	    {
1267	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000040f3 <+16627>:	mulps  %xmm0,%xmm12
   0x0000000000004115 <+16661>:	mulps  %xmm0,%xmm11
   0x000000000000412d <+16685>:	mulps  %xmm0,%xmm10
   0x0000000000004145 <+16709>:	mulps  %xmm0,%xmm9

1268	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_33]), C_row[i]);
   0x00000000000040f7 <+16631>:	addps  0x280(%r9),%xmm12
   0x0000000000004119 <+16665>:	addps  0x290(%r9),%xmm11
   0x0000000000004131 <+16689>:	addps  0x2a0(%r9),%xmm10
   0x0000000000004149 <+16713>:	addps  0x2b0(%r9),%xmm9

1269	      _mm_store_ps(&C[i*4+C_OFFSET_33], C_row[i]);
   0x0000000000004109 <+16649>:	movaps %xmm12,0x280(%r9)
   0x0000000000004121 <+16673>:	movaps %xmm11,0x290(%r9)
   0x0000000000004139 <+16697>:	movaps %xmm10,0x2a0(%r9)
   0x0000000000004151 <+16721>:	movaps %xmm9,0x2b0(%r9)

1270	    }
1271
1272	    /* Reset C(3,4) matrix accumulators */
1273	    C_row[0] = _mm_setzero_ps();
   0x0000000000004111 <+16657>:	xorps  %xmm12,%xmm12

1274	    C_row[1] = _mm_setzero_ps();
   0x0000000000004129 <+16681>:	xorps  %xmm11,%xmm11

1275	    C_row[2] = _mm_setzero_ps();
   0x0000000000004141 <+16705>:	xorps  %xmm10,%xmm10

1276	    C_row[3] = _mm_setzero_ps();
   0x0000000000004159 <+16729>:	xorps  %xmm9,%xmm9

1277
1278	    if (norm_product[8][3] &&
   0x00000000000040ff <+16639>:	movss  0x188(%rsp),%xmm13
   0x000000000000415d <+16733>:	comiss %xmm1,%xmm13
   0x0000000000004161 <+16737>:	jb     0x468a <stream_kernel+18058>

1279	        norm_product[9][7] &&
   0x0000000000004167 <+16743>:	movss  0x68(%rsp),%xmm13
   0x000000000000416e <+16750>:	comiss %xmm1,%xmm13
   0x0000000000004172 <+16754>:	jb     0x468a <stream_kernel+18058>

1280	        norm_product[10][11] &&
   0x0000000000004178 <+16760>:	movss  0x48(%rsp),%xmm13
   0x000000000000417f <+16767>:	comiss %xmm1,%xmm13
   0x0000000000004183 <+16771>:	jb     0x468a <stream_kernel+18058>

1281	        norm_product[11][15])
   0x0000000000004189 <+16777>:	movss  0x38(%rsp),%xmm13
   0x0000000000004190 <+16784>:	comiss %xmm1,%xmm13
   0x0000000000004194 <+16788>:	jb     0x468a <stream_kernel+18058>

1282	    {
1283	      /* A(3,1)*B(1,4) = C(3,4). */
1284	      for (i = 0; i < 4; i++)
1285	      {
1286	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
1287	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
1288	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000419a <+16794>:	movaps 0x800(%r8),%xmm14
   0x00000000000041ba <+16826>:	movaps 0x840(%r8),%xmm10
   0x00000000000041d2 <+16850>:	mulps  0xc0(%rcx),%xmm14
   0x00000000000041f2 <+16882>:	mulps  0xc0(%rcx),%xmm10
   0x000000000000420e <+16910>:	movaps 0x8c0(%r8),%xmm14
   0x0000000000004216 <+16918>:	mulps  0xc0(%rcx),%xmm14
   0x0000000000004236 <+16950>:	movaps 0x880(%r8),%xmm11
   0x000000000000423e <+16958>:	mulps  0xc0(%rcx),%xmm11

1289
1290	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
1291	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
1292	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000041a2 <+16802>:	movaps 0x810(%r8),%xmm12
   0x00000000000041c2 <+16834>:	movaps 0x850(%r8),%xmm9
   0x00000000000041da <+16858>:	mulps  0xd0(%rcx),%xmm12
   0x00000000000041fa <+16890>:	mulps  0xd0(%rcx),%xmm9
   0x000000000000420a <+16906>:	addps  %xmm14,%xmm12
   0x0000000000004246 <+16966>:	addps  %xmm10,%xmm9
   0x000000000000424a <+16970>:	movaps 0x890(%r8),%xmm10
   0x0000000000004252 <+16978>:	mulps  0xd0(%rcx),%xmm10
   0x0000000000004282 <+17026>:	addps  %xmm11,%xmm10
   0x000000000000429a <+17050>:	movaps 0x8d0(%r8),%xmm10
   0x00000000000042a2 <+17058>:	mulps  0xd0(%rcx),%xmm10
   0x00000000000042be <+17086>:	addps  %xmm14,%xmm10

1293
1294	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
1295	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
1296	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000041aa <+16810>:	movaps 0x820(%r8),%xmm11
   0x00000000000041ca <+16842>:	movaps 0x860(%r8),%xmm15
   0x00000000000041e2 <+16866>:	mulps  0xe0(%rcx),%xmm11
   0x0000000000004202 <+16898>:	mulps  0xe0(%rcx),%xmm15
   0x000000000000421e <+16926>:	addps  %xmm12,%xmm11
   0x000000000000425a <+16986>:	addps  %xmm9,%xmm15
   0x000000000000425e <+16990>:	movaps 0x8a0(%r8),%xmm9
   0x0000000000004266 <+16998>:	mulps  0xe0(%rcx),%xmm9
   0x0000000000004296 <+17046>:	addps  %xmm10,%xmm9
   0x00000000000042ae <+17070>:	movaps 0x8e0(%r8),%xmm9
   0x00000000000042b6 <+17078>:	mulps  0xe0(%rcx),%xmm9
   0x00000000000042d2 <+17106>:	addps  %xmm10,%xmm9

1297
1298	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
1299	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
1300	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000041b2 <+16818>:	movaps 0x830(%r8),%xmm13
   0x00000000000041ea <+16874>:	mulps  0xf0(%rcx),%xmm13
   0x0000000000004222 <+16930>:	movaps 0x870(%r8),%xmm12
   0x000000000000422a <+16938>:	mulps  0xf0(%rcx),%xmm12
   0x0000000000004232 <+16946>:	addps  %xmm11,%xmm13
   0x000000000000426e <+17006>:	addps  %xmm15,%xmm12
   0x0000000000004286 <+17030>:	movaps 0x8b0(%r8),%xmm11
   0x000000000000428e <+17038>:	mulps  0xf0(%rcx),%xmm11
   0x00000000000042aa <+17066>:	addps  %xmm9,%xmm11
   0x00000000000042d6 <+17110>:	movaps 0x8f0(%r8),%xmm10
   0x00000000000042de <+17118>:	mulps  0xf0(%rcx),%xmm10
   0x00000000000042fa <+17146>:	addps  %xmm9,%xmm10

1301	      }
1302
1303	      /* A(3,2)*B(2,4) = C(3,4). */
1304	      for (i = 0; i < 4; i++)
1305	      {
1306	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1307	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
1308	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004272 <+17010>:	movaps 0x900(%r8),%xmm15
   0x000000000000427a <+17018>:	mulps  0x1c0(%rcx),%xmm15
   0x00000000000042e6 <+17126>:	addps  %xmm13,%xmm15
   0x0000000000004312 <+17170>:	movaps 0x940(%r8),%xmm15
   0x000000000000431a <+17178>:	mulps  0x1c0(%rcx),%xmm15
   0x000000000000433a <+17210>:	movaps 0x980(%r8),%xmm14
   0x0000000000004342 <+17218>:	mulps  0x1c0(%rcx),%xmm14
   0x000000000000434a <+17226>:	addps  %xmm12,%xmm15
   0x000000000000435e <+17246>:	addps  %xmm11,%xmm14
   0x000000000000439e <+17310>:	movaps 0x9c0(%r8),%xmm13
   0x00000000000043a6 <+17318>:	mulps  0x1c0(%rcx),%xmm13
   0x00000000000043ea <+17386>:	addps  %xmm10,%xmm13

1309
1310	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1311	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
1312	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000042ea <+17130>:	movaps 0x910(%r8),%xmm13
   0x00000000000042f2 <+17138>:	mulps  0x1d0(%rcx),%xmm13
   0x000000000000430e <+17166>:	addps  %xmm15,%xmm13
   0x000000000000434e <+17230>:	movaps 0x950(%r8),%xmm12
   0x0000000000004356 <+17238>:	mulps  0x1d0(%rcx),%xmm12
   0x0000000000004362 <+17250>:	movaps 0x990(%r8),%xmm11
   0x000000000000436a <+17258>:	mulps  0x1d0(%rcx),%xmm11
   0x0000000000004372 <+17266>:	addps  %xmm15,%xmm12
   0x00000000000043ae <+17326>:	addps  %xmm14,%xmm11
   0x00000000000043ee <+17390>:	movaps 0x9d0(%r8),%xmm10
   0x00000000000043f6 <+17398>:	mulps  0x1d0(%rcx),%xmm10
   0x0000000000004412 <+17426>:	addps  %xmm13,%xmm10

1313
1314	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1315	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
1316	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000042c2 <+17090>:	movaps 0x920(%r8),%xmm14
   0x00000000000042ca <+17098>:	mulps  0x1e0(%rcx),%xmm14
   0x0000000000004322 <+17186>:	addps  %xmm13,%xmm14
   0x0000000000004326 <+17190>:	movaps 0x960(%r8),%xmm13
   0x000000000000432e <+17198>:	mulps  0x1e0(%rcx),%xmm13
   0x0000000000004386 <+17286>:	addps  %xmm12,%xmm13
   0x000000000000438a <+17290>:	movaps 0x9a0(%r8),%xmm12
   0x0000000000004392 <+17298>:	mulps  0x1e0(%rcx),%xmm12
   0x00000000000043c2 <+17346>:	addps  %xmm11,%xmm12
   0x00000000000043da <+17370>:	movaps 0x9e0(%r8),%xmm12
   0x00000000000043e2 <+17378>:	mulps  0x1e0(%rcx),%xmm12
   0x0000000000004426 <+17446>:	addps  %xmm10,%xmm12

1317
1318	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1319	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
1320	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000042fe <+17150>:	movaps 0x930(%r8),%xmm9
   0x0000000000004306 <+17158>:	mulps  0x1f0(%rcx),%xmm9
   0x0000000000004336 <+17206>:	addps  %xmm14,%xmm9
   0x0000000000004376 <+17270>:	movaps 0x970(%r8),%xmm15
   0x000000000000437e <+17278>:	mulps  0x1f0(%rcx),%xmm15
   0x000000000000439a <+17306>:	addps  %xmm13,%xmm15
   0x00000000000043b2 <+17330>:	movaps 0x9b0(%r8),%xmm14
   0x00000000000043ba <+17338>:	mulps  0x1f0(%rcx),%xmm14
   0x00000000000043d6 <+17366>:	addps  %xmm12,%xmm14
   0x0000000000004416 <+17430>:	movaps 0x9f0(%r8),%xmm13
   0x000000000000441e <+17438>:	mulps  0x1f0(%rcx),%xmm13
   0x000000000000443a <+17466>:	addps  %xmm12,%xmm13

1321	      }
1322
1323	      /* A(3,3)*B(3,4) = C(3,4). */
1324	      for (i = 0; i < 4; i++)
1325	      {
1326	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1327	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
1328	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000043c6 <+17350>:	movaps 0xa00(%r8),%xmm11
   0x00000000000043ce <+17358>:	mulps  0x2c0(%rcx),%xmm11
   0x00000000000043fe <+17406>:	addps  %xmm9,%xmm11
   0x000000000000442a <+17450>:	movaps 0xa40(%r8),%xmm10
   0x0000000000004432 <+17458>:	mulps  0x2c0(%rcx),%xmm10
   0x0000000000004462 <+17506>:	addps  %xmm15,%xmm10
   0x000000000000447a <+17530>:	movaps 0xa80(%r8),%xmm9
   0x0000000000004482 <+17538>:	mulps  0x2c0(%rcx),%xmm9
   0x00000000000044b2 <+17586>:	addps  %xmm14,%xmm9
   0x00000000000044ca <+17610>:	movaps 0xac0(%r8),%xmm15
   0x00000000000044d2 <+17618>:	mulps  0x2c0(%rcx),%xmm15
   0x0000000000004502 <+17666>:	addps  %xmm13,%xmm15

1329
1330	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1331	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
1332	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004402 <+17410>:	movaps 0xa10(%r8),%xmm9
   0x000000000000440a <+17418>:	mulps  0x2d0(%rcx),%xmm9
   0x000000000000444e <+17486>:	addps  %xmm11,%xmm9
   0x0000000000004466 <+17510>:	movaps 0xa50(%r8),%xmm15
   0x000000000000446e <+17518>:	mulps  0x2d0(%rcx),%xmm15
   0x000000000000449e <+17566>:	addps  %xmm10,%xmm15
   0x00000000000044b6 <+17590>:	movaps 0xa90(%r8),%xmm14
   0x00000000000044be <+17598>:	mulps  0x2d0(%rcx),%xmm14
   0x00000000000044ee <+17646>:	addps  %xmm9,%xmm14
   0x0000000000004506 <+17670>:	movaps 0xad0(%r8),%xmm13
   0x000000000000450e <+17678>:	mulps  0x2d0(%rcx),%xmm13
   0x000000000000453e <+17726>:	addps  %xmm15,%xmm13

1333
1334	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1335	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
1336	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004452 <+17490>:	movaps 0xa20(%r8),%xmm11
   0x000000000000445a <+17498>:	mulps  0x2e0(%rcx),%xmm11
   0x0000000000004476 <+17526>:	addps  %xmm9,%xmm11
   0x00000000000044a2 <+17570>:	movaps 0xa60(%r8),%xmm10
   0x00000000000044aa <+17578>:	mulps  0x2e0(%rcx),%xmm10
   0x00000000000044c6 <+17606>:	addps  %xmm15,%xmm10
   0x00000000000044f2 <+17650>:	movaps 0xaa0(%r8),%xmm9
   0x00000000000044fa <+17658>:	mulps  0x2e0(%rcx),%xmm9
   0x0000000000004516 <+17686>:	addps  %xmm14,%xmm9
   0x000000000000451a <+17690>:	movaps 0xae0(%r8),%xmm14
   0x0000000000004522 <+17698>:	mulps  0x2e0(%rcx),%xmm14
   0x0000000000004552 <+17746>:	addps  %xmm13,%xmm14

1337
1338	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1339	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
1340	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000443e <+17470>:	movaps 0xa30(%r8),%xmm12
   0x0000000000004446 <+17478>:	mulps  0x2f0(%rcx),%xmm12
   0x000000000000448a <+17546>:	addps  %xmm11,%xmm12
   0x000000000000448e <+17550>:	movaps 0xa70(%r8),%xmm11
   0x0000000000004496 <+17558>:	mulps  0x2f0(%rcx),%xmm11
   0x00000000000044da <+17626>:	addps  %xmm10,%xmm11
   0x00000000000044de <+17630>:	movaps 0xab0(%r8),%xmm10
   0x00000000000044e6 <+17638>:	mulps  0x2f0(%rcx),%xmm10
   0x000000000000452a <+17706>:	addps  %xmm9,%xmm10
   0x000000000000452e <+17710>:	movaps 0xaf0(%r8),%xmm9
   0x0000000000004536 <+17718>:	mulps  0x2f0(%rcx),%xmm9
   0x0000000000004566 <+17766>:	addps  %xmm14,%xmm9

1341	      }
1342
1343	      /* A(3,4)*B(4,4) = C(3,4). */
1344	      for (i = 0; i < 4; i++)
1345	      {
1346	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1347	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
1348	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004542 <+17730>:	movaps 0xb00(%r8),%xmm15
   0x000000000000454a <+17738>:	mulps  0x3c0(%rcx),%xmm15
   0x000000000000456a <+17770>:	movaps 0xb40(%r8),%xmm14
   0x0000000000004572 <+17778>:	mulps  0x3c0(%rcx),%xmm14
   0x000000000000457a <+17786>:	addps  %xmm12,%xmm15
   0x000000000000458e <+17806>:	addps  %xmm11,%xmm14
   0x00000000000045a6 <+17830>:	movaps 0xb80(%r8),%xmm15
   0x00000000000045ae <+17838>:	mulps  0x3c0(%rcx),%xmm15
   0x00000000000045ce <+17870>:	movaps 0xbc0(%r8),%xmm14
   0x00000000000045d6 <+17878>:	mulps  0x3c0(%rcx),%xmm14
   0x00000000000045f2 <+17906>:	addps  %xmm10,%xmm15
   0x000000000000461e <+17950>:	addps  %xmm9,%xmm14

1349
1350	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1351	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
1352	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000457e <+17790>:	movaps 0xb10(%r8),%xmm12
   0x0000000000004586 <+17798>:	mulps  0x3d0(%rcx),%xmm12
   0x0000000000004592 <+17810>:	movaps 0xb50(%r8),%xmm11
   0x000000000000459a <+17818>:	mulps  0x3d0(%rcx),%xmm11
   0x00000000000045a2 <+17826>:	addps  %xmm15,%xmm12
   0x00000000000045ca <+17866>:	addps  %xmm14,%xmm11
   0x00000000000045f6 <+17910>:	movaps 0xb90(%r8),%xmm10
   0x00000000000045fe <+17918>:	mulps  0x3d0(%rcx),%xmm10
   0x000000000000461a <+17946>:	addps  %xmm15,%xmm10
   0x0000000000004622 <+17954>:	movaps 0xbd0(%r8),%xmm9
   0x000000000000462a <+17962>:	mulps  0x3d0(%rcx),%xmm9
   0x0000000000004646 <+17990>:	addps  %xmm14,%xmm9

1353
1354	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1355	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
1356	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004556 <+17750>:	movaps 0xb20(%r8),%xmm13
   0x000000000000455e <+17758>:	mulps  0x3e0(%rcx),%xmm13
   0x00000000000045b6 <+17846>:	addps  %xmm12,%xmm13
   0x00000000000045e2 <+17890>:	movaps 0xb60(%r8),%xmm13
   0x00000000000045ea <+17898>:	mulps  0x3e0(%rcx),%xmm13
   0x0000000000004606 <+17926>:	addps  %xmm11,%xmm13
   0x0000000000004636 <+17974>:	movaps 0xba0(%r8),%xmm13
   0x000000000000463e <+17982>:	mulps  0x3e0(%rcx),%xmm13
   0x000000000000464a <+17994>:	addps  %xmm10,%xmm13
   0x0000000000004662 <+18018>:	movaps 0xbe0(%r8),%xmm13
   0x000000000000466a <+18026>:	mulps  0x3e0(%rcx),%xmm13
   0x0000000000004672 <+18034>:	addps  %xmm9,%xmm13

1357
1358	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1359	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
1360	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000045ba <+17850>:	movaps 0xb30(%r8),%xmm12
   0x00000000000045c2 <+17858>:	mulps  0x3f0(%rcx),%xmm12
   0x00000000000045de <+17886>:	addps  %xmm13,%xmm12
   0x000000000000460a <+17930>:	movaps 0xb70(%r8),%xmm11
   0x0000000000004612 <+17938>:	mulps  0x3f0(%rcx),%xmm11
   0x0000000000004632 <+17970>:	addps  %xmm13,%xmm11
   0x000000000000464e <+17998>:	movaps 0xbb0(%r8),%xmm10
   0x0000000000004656 <+18006>:	mulps  0x3f0(%rcx),%xmm10
   0x000000000000465e <+18014>:	addps  %xmm13,%xmm10
   0x0000000000004676 <+18038>:	movaps 0xbf0(%r8),%xmm9
   0x000000000000467e <+18046>:	mulps  0x3f0(%rcx),%xmm9
   0x0000000000004686 <+18054>:	addps  %xmm13,%xmm9

1361	      }
1362	    }
1363
1364	    /* Store C(3,4) block. */
1365	    for (i = 0; i < 4; i++)
1366	    {
1367	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x000000000000468a <+18058>:	mulps  %xmm0,%xmm12
   0x00000000000046a2 <+18082>:	mulps  %xmm0,%xmm11
   0x00000000000046ba <+18106>:	mulps  %xmm0,%xmm10
   0x00000000000046d2 <+18130>:	mulps  %xmm0,%xmm9

1368	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_34]), C_row[i]);
   0x000000000000468e <+18062>:	addps  0x2c0(%r9),%xmm12
   0x00000000000046a6 <+18086>:	addps  0x2d0(%r9),%xmm11
   0x00000000000046be <+18110>:	addps  0x2e0(%r9),%xmm10
   0x00000000000046d6 <+18134>:	addps  0x2f0(%r9),%xmm9

1369	      _mm_store_ps(&C[i*4+C_OFFSET_34], C_row[i]);
   0x0000000000004696 <+18070>:	movaps %xmm12,0x2c0(%r9)
   0x00000000000046ae <+18094>:	movaps %xmm11,0x2d0(%r9)
   0x00000000000046c6 <+18118>:	movaps %xmm10,0x2e0(%r9)
   0x00000000000046de <+18142>:	movaps %xmm9,0x2f0(%r9)

1370	    }
1371
1372	    /* Reset C(4,1) matrix accumulators */
1373	    C_row[0] = _mm_setzero_ps();
   0x000000000000469e <+18078>:	xorps  %xmm12,%xmm12

1374	    C_row[1] = _mm_setzero_ps();
   0x00000000000046b6 <+18102>:	xorps  %xmm11,%xmm11

1375	    C_row[2] = _mm_setzero_ps();
   0x00000000000046ce <+18126>:	xorps  %xmm10,%xmm10

1376	    C_row[3] = _mm_setzero_ps();
   0x00000000000046e6 <+18150>:	xorps  %xmm9,%xmm9

1377
1378	    if (norm_product[12][0] &&
   0x00000000000046ea <+18154>:	comiss %xmm1,%xmm8
   0x00000000000046ee <+18158>:	jb     0x4be6 <stream_kernel+19430>

1379	        norm_product[13][4] &&
   0x00000000000046f4 <+18164>:	movss  0x90(%rsp),%xmm8
   0x00000000000046fe <+18174>:	comiss %xmm1,%xmm8
   0x0000000000004702 <+18178>:	jb     0x4be6 <stream_kernel+19430>

1380	        norm_product[14][8] &&
   0x0000000000004708 <+18184>:	movss  0x40(%rsp),%xmm8
   0x000000000000470f <+18191>:	comiss %xmm1,%xmm8
   0x0000000000004713 <+18195>:	jb     0x4be6 <stream_kernel+19430>

1381	        norm_product[15][12])
   0x0000000000004719 <+18201>:	movss  0x28(%rsp),%xmm8
   0x0000000000004720 <+18208>:	comiss %xmm1,%xmm8
   0x0000000000004724 <+18212>:	jb     0x4be6 <stream_kernel+19430>

1382	    {
1383	      /* A(4,1)*B(1,1) = C(4,1). */
1384	      for (i = 0; i < 4; i++)
1385	      {
1386	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1387	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
1388	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000472a <+18218>:	movaps 0xc00(%r8),%xmm12
   0x000000000000475a <+18266>:	movaps 0xc80(%r8),%xmm8
   0x000000000000476a <+18282>:	mulps  (%rcx),%xmm12
   0x0000000000004787 <+18311>:	mulps  (%rcx),%xmm8
   0x0000000000004794 <+18324>:	movaps 0xc40(%r8),%xmm12
   0x000000000000479c <+18332>:	mulps  (%rcx),%xmm12
   0x00000000000047c9 <+18377>:	movaps 0xcc0(%r8),%xmm12
   0x00000000000047d1 <+18385>:	mulps  (%rcx),%xmm12

1389
1390	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1391	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
1392	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004732 <+18226>:	movaps 0xc10(%r8),%xmm13
   0x000000000000474a <+18250>:	movaps 0xc50(%r8),%xmm10
   0x000000000000476e <+18286>:	mulps  0x10(%rcx),%xmm13
   0x000000000000477d <+18301>:	mulps  0x10(%rcx),%xmm10
   0x0000000000004790 <+18320>:	addps  %xmm12,%xmm13
   0x00000000000047a4 <+18340>:	movaps 0xc90(%r8),%xmm13
   0x00000000000047ac <+18348>:	mulps  0x10(%rcx),%xmm13
   0x00000000000047c5 <+18373>:	addps  %xmm12,%xmm10
   0x00000000000047e6 <+18406>:	addps  %xmm8,%xmm13
   0x00000000000047fb <+18427>:	movaps 0xcd0(%r8),%xmm9
   0x0000000000004803 <+18435>:	mulps  0x10(%rcx),%xmm9
   0x000000000000482d <+18477>:	addps  %xmm12,%xmm9

1393
1394	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1395	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
1396	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000473a <+18234>:	movaps 0xc20(%r8),%xmm14
   0x0000000000004752 <+18258>:	movaps 0xc60(%r8),%xmm9
   0x0000000000004762 <+18274>:	movaps 0xca0(%r8),%xmm15
   0x0000000000004773 <+18291>:	mulps  0x20(%rcx),%xmm14
   0x0000000000004782 <+18306>:	mulps  0x20(%rcx),%xmm9
   0x000000000000478b <+18315>:	mulps  0x20(%rcx),%xmm15
   0x00000000000047a0 <+18336>:	addps  %xmm13,%xmm14
   0x00000000000047d5 <+18389>:	addps  %xmm10,%xmm9
   0x0000000000004808 <+18440>:	addps  %xmm13,%xmm15
   0x000000000000480c <+18444>:	movaps 0xce0(%r8),%xmm13
   0x0000000000004814 <+18452>:	mulps  0x20(%rcx),%xmm13
   0x0000000000004841 <+18497>:	addps  %xmm9,%xmm13

1397
1398	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1399	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
1400	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004742 <+18242>:	movaps 0xc30(%r8),%xmm11
   0x0000000000004778 <+18296>:	mulps  0x30(%rcx),%xmm11
   0x00000000000047b1 <+18353>:	addps  %xmm14,%xmm11
   0x00000000000047d9 <+18393>:	movaps 0xc70(%r8),%xmm10
   0x00000000000047e1 <+18401>:	mulps  0x30(%rcx),%xmm10
   0x00000000000047ea <+18410>:	movaps 0xcb0(%r8),%xmm8
   0x00000000000047f2 <+18418>:	mulps  0x30(%rcx),%xmm8
   0x00000000000047f7 <+18423>:	addps  %xmm9,%xmm10
   0x0000000000004819 <+18457>:	addps  %xmm15,%xmm8
   0x0000000000004845 <+18501>:	movaps 0xcf0(%r8),%xmm9
   0x000000000000484d <+18509>:	mulps  0x30(%rcx),%xmm9
   0x0000000000004866 <+18534>:	addps  %xmm13,%xmm9

1401	      }
1402
1403	      /* A(4,2)*B(2,1) = C(4,1). */
1404	      for (i = 0; i < 4; i++)
1405	      {
1406	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1407	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
1408	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000047b5 <+18357>:	movaps 0xd00(%r8),%xmm14
   0x00000000000047bd <+18365>:	mulps  0x100(%rcx),%xmm14
   0x0000000000004852 <+18514>:	addps  %xmm11,%xmm14
   0x000000000000486a <+18538>:	movaps 0xd40(%r8),%xmm13
   0x0000000000004872 <+18546>:	mulps  0x100(%rcx),%xmm13
   0x0000000000004892 <+18578>:	movaps 0xd80(%r8),%xmm11
   0x000000000000489a <+18586>:	mulps  0x100(%rcx),%xmm11
   0x00000000000048b6 <+18614>:	addps  %xmm10,%xmm13
   0x00000000000048ca <+18634>:	addps  %xmm8,%xmm11
   0x000000000000491e <+18718>:	movaps 0xdc0(%r8),%xmm11
   0x0000000000004926 <+18726>:	mulps  0x100(%rcx),%xmm11
   0x0000000000004956 <+18774>:	addps  %xmm9,%xmm11

1409
1410	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1411	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
1412	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004856 <+18518>:	movaps 0xd10(%r8),%xmm11
   0x000000000000485e <+18526>:	mulps  0x110(%rcx),%xmm11
   0x000000000000487a <+18554>:	addps  %xmm14,%xmm11
   0x00000000000048ba <+18618>:	movaps 0xd50(%r8),%xmm10
   0x00000000000048c2 <+18626>:	mulps  0x110(%rcx),%xmm10
   0x00000000000048ce <+18638>:	movaps 0xd90(%r8),%xmm8
   0x00000000000048d6 <+18646>:	mulps  0x110(%rcx),%xmm8
   0x00000000000048de <+18654>:	addps  %xmm13,%xmm10
   0x000000000000491a <+18714>:	addps  %xmm11,%xmm8
   0x000000000000495a <+18778>:	movaps 0xdd0(%r8),%xmm9
   0x0000000000004962 <+18786>:	mulps  0x110(%rcx),%xmm9
   0x000000000000497e <+18814>:	addps  %xmm11,%xmm9

1413
1414	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1415	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
1416	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004831 <+18481>:	movaps 0xd20(%r8),%xmm12
   0x0000000000004839 <+18489>:	mulps  0x120(%rcx),%xmm12
   0x000000000000488e <+18574>:	addps  %xmm11,%xmm12
   0x00000000000048a6 <+18598>:	movaps 0xd60(%r8),%xmm12
   0x00000000000048ae <+18606>:	mulps  0x120(%rcx),%xmm12
   0x00000000000048f2 <+18674>:	addps  %xmm10,%xmm12
   0x000000000000490a <+18698>:	movaps 0xda0(%r8),%xmm12
   0x0000000000004912 <+18706>:	mulps  0x120(%rcx),%xmm12
   0x000000000000492e <+18734>:	addps  %xmm8,%xmm12
   0x0000000000004982 <+18818>:	movaps 0xde0(%r8),%xmm11
   0x000000000000498a <+18826>:	mulps  0x120(%rcx),%xmm11
   0x00000000000049a6 <+18854>:	addps  %xmm9,%xmm11

1417
1418	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1419	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
1420	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000481d <+18461>:	movaps 0xd30(%r8),%xmm15
   0x0000000000004825 <+18469>:	mulps  0x130(%rcx),%xmm15
   0x000000000000487e <+18558>:	movaps 0xd70(%r8),%xmm14
   0x0000000000004886 <+18566>:	mulps  0x130(%rcx),%xmm14
   0x00000000000048a2 <+18594>:	addps  %xmm12,%xmm15
   0x00000000000048e2 <+18658>:	movaps 0xdb0(%r8),%xmm13
   0x00000000000048ea <+18666>:	mulps  0x130(%rcx),%xmm13
   0x0000000000004906 <+18694>:	addps  %xmm12,%xmm14
   0x0000000000004942 <+18754>:	addps  %xmm12,%xmm13
   0x0000000000004946 <+18758>:	movaps 0xdf0(%r8),%xmm12
   0x000000000000494e <+18766>:	mulps  0x130(%rcx),%xmm12
   0x00000000000049ba <+18874>:	addps  %xmm11,%xmm12

1421	      }
1422
1423	      /* A(4,3)*B(3,1) = C(4,1). */
1424	      for (i = 0; i < 4; i++)
1425	      {
1426	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1427	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
1428	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000048f6 <+18678>:	movaps 0xe00(%r8),%xmm10
   0x00000000000048fe <+18686>:	mulps  0x200(%rcx),%xmm10
   0x0000000000004932 <+18738>:	movaps 0xe80(%r8),%xmm8
   0x000000000000493a <+18746>:	mulps  0x200(%rcx),%xmm8
   0x000000000000496a <+18794>:	addps  %xmm15,%xmm10
   0x00000000000049aa <+18858>:	movaps 0xe40(%r8),%xmm9
   0x00000000000049b2 <+18866>:	mulps  0x200(%rcx),%xmm9
   0x00000000000049f6 <+18934>:	addps  %xmm14,%xmm9
   0x0000000000004a0a <+18954>:	addps  %xmm13,%xmm8
   0x0000000000004a4a <+19018>:	movaps 0xec0(%r8),%xmm14
   0x0000000000004a52 <+19026>:	mulps  0x200(%rcx),%xmm14
   0x0000000000004a6e <+19054>:	addps  %xmm12,%xmm14

1429
1430	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1431	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
1432	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000496e <+18798>:	movaps 0xe10(%r8),%xmm15
   0x0000000000004976 <+18806>:	mulps  0x210(%rcx),%xmm15
   0x0000000000004992 <+18834>:	addps  %xmm10,%xmm15
   0x00000000000049fa <+18938>:	movaps 0xe50(%r8),%xmm14
   0x0000000000004a02 <+18946>:	mulps  0x210(%rcx),%xmm14
   0x0000000000004a0e <+18958>:	movaps 0xe90(%r8),%xmm13
   0x0000000000004a16 <+18966>:	mulps  0x210(%rcx),%xmm13
   0x0000000000004a1e <+18974>:	addps  %xmm9,%xmm14
   0x0000000000004a32 <+18994>:	addps  %xmm8,%xmm13
   0x0000000000004a72 <+19058>:	movaps 0xed0(%r8),%xmm12
   0x0000000000004a7a <+19066>:	mulps  0x210(%rcx),%xmm12
   0x0000000000004aaa <+19114>:	addps  %xmm14,%xmm12

1433
1434	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1435	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
1436	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004996 <+18838>:	movaps 0xe20(%r8),%xmm10
   0x000000000000499e <+18846>:	mulps  0x220(%rcx),%xmm10
   0x00000000000049ce <+18894>:	addps  %xmm15,%xmm10
   0x0000000000004a22 <+18978>:	movaps 0xe60(%r8),%xmm9
   0x0000000000004a2a <+18986>:	mulps  0x220(%rcx),%xmm9
   0x0000000000004a46 <+19014>:	addps  %xmm14,%xmm9
   0x0000000000004a5e <+19038>:	movaps 0xea0(%r8),%xmm9
   0x0000000000004a66 <+19046>:	mulps  0x220(%rcx),%xmm9
   0x0000000000004a82 <+19074>:	addps  %xmm13,%xmm9
   0x0000000000004a86 <+19078>:	movaps 0xee0(%r8),%xmm13
   0x0000000000004a8e <+19086>:	mulps  0x220(%rcx),%xmm13
   0x0000000000004abe <+19134>:	addps  %xmm12,%xmm13

1437
1438	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1439	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
1440	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000049be <+18878>:	movaps 0xe30(%r8),%xmm11
   0x00000000000049c6 <+18886>:	mulps  0x230(%rcx),%xmm11
   0x00000000000049e2 <+18914>:	addps  %xmm10,%xmm11
   0x00000000000049e6 <+18918>:	movaps 0xe70(%r8),%xmm10
   0x00000000000049ee <+18926>:	mulps  0x230(%rcx),%xmm10
   0x0000000000004a36 <+18998>:	movaps 0xeb0(%r8),%xmm8
   0x0000000000004a3e <+19006>:	mulps  0x230(%rcx),%xmm8
   0x0000000000004a5a <+19034>:	addps  %xmm9,%xmm10
   0x0000000000004a96 <+19094>:	addps  %xmm9,%xmm8
   0x0000000000004a9a <+19098>:	movaps 0xef0(%r8),%xmm9
   0x0000000000004aa2 <+19106>:	mulps  0x230(%rcx),%xmm9
   0x0000000000004ad2 <+19154>:	addps  %xmm13,%xmm9

1441	      }
1442
1443	      /* A(4,4)*B(4,1) = C(4,1). */
1444	      for (i = 0; i < 4; i++)
1445	      {
1446	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1447	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
1448	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004aae <+19118>:	movaps 0xf00(%r8),%xmm14
   0x0000000000004ab6 <+19126>:	mulps  0x300(%rcx),%xmm14
   0x0000000000004ae6 <+19174>:	addps  %xmm11,%xmm14
   0x0000000000004afe <+19198>:	movaps 0xf80(%r8),%xmm14
   0x0000000000004b06 <+19206>:	mulps  0x300(%rcx),%xmm14
   0x0000000000004b12 <+19218>:	movaps 0xf40(%r8),%xmm11
   0x0000000000004b1a <+19226>:	mulps  0x300(%rcx),%xmm11
   0x0000000000004b36 <+19254>:	addps  %xmm10,%xmm11
   0x0000000000004b4a <+19274>:	addps  %xmm8,%xmm14
   0x0000000000004b9e <+19358>:	movaps 0xfc0(%r8),%xmm14
   0x0000000000004ba6 <+19366>:	mulps  0x300(%rcx),%xmm14
   0x0000000000004bb6 <+19382>:	addps  %xmm9,%xmm14

1449
1450	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1451	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
1452	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004aea <+19178>:	movaps 0xf10(%r8),%xmm11
   0x0000000000004af2 <+19186>:	mulps  0x310(%rcx),%xmm11
   0x0000000000004afa <+19194>:	addps  %xmm14,%xmm11
   0x0000000000004b3a <+19258>:	movaps 0xf50(%r8),%xmm10
   0x0000000000004b42 <+19266>:	mulps  0x310(%rcx),%xmm10
   0x0000000000004b4e <+19278>:	movaps 0xf90(%r8),%xmm8
   0x0000000000004b56 <+19286>:	mulps  0x310(%rcx),%xmm8
   0x0000000000004b5e <+19294>:	addps  %xmm11,%xmm10
   0x0000000000004b9a <+19354>:	addps  %xmm14,%xmm8
   0x0000000000004bba <+19386>:	movaps 0xfd0(%r8),%xmm9
   0x0000000000004bc2 <+19394>:	mulps  0x310(%rcx),%xmm9
   0x0000000000004bca <+19402>:	addps  %xmm14,%xmm9

1453
1454	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1455	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
1456	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000049d2 <+18898>:	movaps 0xf20(%r8),%xmm15
   0x00000000000049da <+18906>:	mulps  0x320(%rcx),%xmm15
   0x0000000000004ad6 <+19158>:	movaps 0xf60(%r8),%xmm13
   0x0000000000004ade <+19166>:	mulps  0x320(%rcx),%xmm13
   0x0000000000004b0e <+19214>:	addps  %xmm11,%xmm15
   0x0000000000004b26 <+19238>:	movaps 0xfe0(%r8),%xmm15
   0x0000000000004b2e <+19246>:	mulps  0x320(%rcx),%xmm15
   0x0000000000004b72 <+19314>:	addps  %xmm10,%xmm13
   0x0000000000004b8a <+19338>:	movaps 0xfa0(%r8),%xmm13
   0x0000000000004b92 <+19346>:	mulps  0x320(%rcx),%xmm13
   0x0000000000004bae <+19374>:	addps  %xmm8,%xmm13
   0x0000000000004bce <+19406>:	addps  %xmm9,%xmm15

1457
1458	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1459	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
1460	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004ac2 <+19138>:	movaps 0xf30(%r8),%xmm12
   0x0000000000004aca <+19146>:	mulps  0x330(%rcx),%xmm12
   0x0000000000004b22 <+19234>:	addps  %xmm15,%xmm12
   0x0000000000004b62 <+19298>:	movaps 0xf70(%r8),%xmm11
   0x0000000000004b6a <+19306>:	mulps  0x330(%rcx),%xmm11
   0x0000000000004b76 <+19318>:	movaps 0xfb0(%r8),%xmm10
   0x0000000000004b7e <+19326>:	mulps  0x330(%rcx),%xmm10
   0x0000000000004b86 <+19334>:	addps  %xmm13,%xmm11
   0x0000000000004bb2 <+19378>:	addps  %xmm13,%xmm10
   0x0000000000004bd2 <+19410>:	movaps 0xff0(%r8),%xmm9
   0x0000000000004bda <+19418>:	mulps  0x330(%rcx),%xmm9
   0x0000000000004be2 <+19426>:	addps  %xmm15,%xmm9

1461	      }
1462	    }
1463
1464	    /* Store C(4,1) block. */
1465	    for (i = 0; i < 4; i++)
1466	    {
1467	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000004be6 <+19430>:	mulps  %xmm0,%xmm12
   0x0000000000004c02 <+19458>:	mulps  %xmm0,%xmm11
   0x0000000000004c16 <+19478>:	mulps  %xmm0,%xmm10
   0x0000000000004c2e <+19502>:	mulps  %xmm0,%xmm9

1468	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_41]), C_row[i]);
   0x0000000000004bea <+19434>:	addps  0x300(%r9),%xmm12
   0x0000000000004c06 <+19462>:	addps  0x310(%r9),%xmm11
   0x0000000000004c1a <+19482>:	addps  0x320(%r9),%xmm10
   0x0000000000004c32 <+19506>:	addps  0x330(%r9),%xmm9

1469	      _mm_store_ps(&C[i*4+C_OFFSET_41], C_row[i]);
   0x0000000000004bf6 <+19446>:	movaps %xmm12,0x300(%r9)
   0x0000000000004c0e <+19470>:	movaps %xmm11,0x310(%r9)
   0x0000000000004c22 <+19490>:	movaps %xmm10,0x320(%r9)
   0x0000000000004c3a <+19514>:	movaps %xmm9,0x330(%r9)

1470	    }
1471
1472	    /* Reset C(4,2) matrix accumulators */
1473	    C_row[0] = _mm_setzero_ps();
   0x0000000000004c2a <+19498>:	xorps  %xmm10,%xmm10

1474	    C_row[1] = _mm_setzero_ps();
   0x0000000000004c42 <+19522>:	xorps  %xmm9,%xmm9

1475	    C_row[2] = _mm_setzero_ps();
   0x0000000000004bf2 <+19442>:	xorps  %xmm8,%xmm8

1476	    C_row[3] = _mm_setzero_ps();
   0x0000000000004bfe <+19454>:	xorps  %xmm12,%xmm12

1477
1478	    if (norm_product[12][1] &&
   0x0000000000004c46 <+19526>:	comiss %xmm1,%xmm7
   0x0000000000004c49 <+19529>:	jb     0x5135 <stream_kernel+20789>

1479	        norm_product[13][5] &&
   0x0000000000004c4f <+19535>:	movss  0x70(%rsp),%xmm7
   0x0000000000004c55 <+19541>:	comiss %xmm1,%xmm7
   0x0000000000004c58 <+19544>:	jb     0x5135 <stream_kernel+20789>

1480	        norm_product[14][9] &&
   0x0000000000004c5e <+19550>:	movss  0x30(%rsp),%xmm7
   0x0000000000004c64 <+19556>:	comiss %xmm1,%xmm7
   0x0000000000004c67 <+19559>:	jb     0x5135 <stream_kernel+20789>

1481	        norm_product[15][13])
   0x0000000000004c6d <+19565>:	movss  0x18(%rsp),%xmm7
   0x0000000000004c73 <+19571>:	comiss %xmm1,%xmm7
   0x0000000000004c76 <+19574>:	jb     0x5135 <stream_kernel+20789>

1482	    {
1483	      /* A(4,1)*B(1,2) = C(4,2). */
1484	      for (i = 0; i < 4; i++)
1485	      {
1486	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1487	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
1488	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004c7c <+19580>:	movaps 0xc00(%r8),%xmm10
   0x0000000000004c9c <+19612>:	movaps 0xc40(%r8),%xmm13
   0x0000000000004cb4 <+19636>:	movaps 0xc80(%r8),%xmm9
   0x0000000000004cbc <+19644>:	mulps  0x40(%rcx),%xmm10
   0x0000000000004ccf <+19663>:	mulps  0x40(%rcx),%xmm13
   0x0000000000004cde <+19678>:	mulps  0x40(%rcx),%xmm9
   0x0000000000004d27 <+19751>:	movaps 0xcc0(%r8),%xmm13
   0x0000000000004d2f <+19759>:	mulps  0x40(%rcx),%xmm13

1489
1490	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1491	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
1492	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004c84 <+19588>:	movaps 0xc10(%r8),%xmm11
   0x0000000000004ca4 <+19620>:	movaps 0xc50(%r8),%xmm14
   0x0000000000004cc1 <+19649>:	mulps  0x50(%rcx),%xmm11
   0x0000000000004cd4 <+19668>:	mulps  0x50(%rcx),%xmm14
   0x0000000000004ce3 <+19683>:	movaps 0xcd0(%r8),%xmm15
   0x0000000000004ceb <+19691>:	addps  %xmm10,%xmm11
   0x0000000000004cef <+19695>:	mulps  0x50(%rcx),%xmm15
   0x0000000000004d00 <+19712>:	movaps 0xc90(%r8),%xmm11
   0x0000000000004d0d <+19725>:	mulps  0x50(%rcx),%xmm11
   0x0000000000004d23 <+19747>:	addps  %xmm13,%xmm14
   0x0000000000004d59 <+19801>:	addps  %xmm9,%xmm11
   0x0000000000004d92 <+19858>:	addps  %xmm13,%xmm15

1493
1494	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1495	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
1496	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004c8c <+19596>:	movaps 0xc20(%r8),%xmm12
   0x0000000000004cc6 <+19654>:	mulps  0x60(%rcx),%xmm12
   0x0000000000004cf4 <+19700>:	movaps 0xc60(%r8),%xmm10
   0x0000000000004cfc <+19708>:	addps  %xmm11,%xmm12
   0x0000000000004d08 <+19720>:	mulps  0x60(%rcx),%xmm10
   0x0000000000004d16 <+19734>:	movaps 0xca0(%r8),%xmm12
   0x0000000000004d1e <+19742>:	mulps  0x60(%rcx),%xmm12
   0x0000000000004d34 <+19764>:	addps  %xmm14,%xmm10
   0x0000000000004d4c <+19788>:	movaps 0xce0(%r8),%xmm10
   0x0000000000004d54 <+19796>:	mulps  0x60(%rcx),%xmm10
   0x0000000000004d6a <+19818>:	addps  %xmm11,%xmm12
   0x0000000000004da6 <+19878>:	addps  %xmm15,%xmm10

1497
1498	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1499	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
1500	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004c94 <+19604>:	movaps 0xc30(%r8),%xmm7
   0x0000000000004cac <+19628>:	movaps 0xc70(%r8),%xmm8
   0x0000000000004ccb <+19659>:	mulps  0x70(%rcx),%xmm7
   0x0000000000004cd9 <+19673>:	mulps  0x70(%rcx),%xmm8
   0x0000000000004d12 <+19730>:	addps  %xmm12,%xmm7
   0x0000000000004d48 <+19784>:	addps  %xmm10,%xmm8
   0x0000000000004d5d <+19805>:	movaps 0xcb0(%r8),%xmm9
   0x0000000000004d65 <+19813>:	mulps  0x70(%rcx),%xmm9
   0x0000000000004d7e <+19838>:	addps  %xmm12,%xmm9
   0x0000000000004daa <+19882>:	movaps 0xcf0(%r8),%xmm15
   0x0000000000004db2 <+19890>:	mulps  0x70(%rcx),%xmm15
   0x0000000000004dca <+19914>:	addps  %xmm10,%xmm15

1501	      }
1502
1503	      /* A(4,2)*B(2,2) = C(4,2). */
1504	      for (i = 0; i < 4; i++)
1505	      {
1506	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1507	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
1508	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004d6e <+19822>:	movaps 0xd00(%r8),%xmm11
   0x0000000000004d76 <+19830>:	mulps  0x140(%rcx),%xmm11
   0x0000000000004db7 <+19895>:	addps  %xmm7,%xmm11
   0x0000000000004de2 <+19938>:	movaps 0xd40(%r8),%xmm11
   0x0000000000004dea <+19946>:	mulps  0x140(%rcx),%xmm11
   0x0000000000004df6 <+19958>:	movaps 0xd80(%r8),%xmm7
   0x0000000000004dfe <+19966>:	mulps  0x140(%rcx),%xmm7
   0x0000000000004e19 <+19993>:	addps  %xmm8,%xmm11
   0x0000000000004e2d <+20013>:	addps  %xmm9,%xmm7
   0x0000000000004e6d <+20077>:	movaps 0xdc0(%r8),%xmm10
   0x0000000000004e75 <+20085>:	mulps  0x140(%rcx),%xmm10
   0x0000000000004eb8 <+20152>:	addps  %xmm15,%xmm10

1509
1510	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1511	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
1512	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004dbb <+19899>:	movaps 0xd10(%r8),%xmm7
   0x0000000000004dc3 <+19907>:	mulps  0x150(%rcx),%xmm7
   0x0000000000004dde <+19934>:	addps  %xmm11,%xmm7
   0x0000000000004e1d <+19997>:	movaps 0xd50(%r8),%xmm8
   0x0000000000004e25 <+20005>:	mulps  0x150(%rcx),%xmm8
   0x0000000000004e31 <+20017>:	movaps 0xd90(%r8),%xmm9
   0x0000000000004e39 <+20025>:	mulps  0x150(%rcx),%xmm9
   0x0000000000004e41 <+20033>:	addps  %xmm11,%xmm8
   0x0000000000004e7d <+20093>:	addps  %xmm7,%xmm9
   0x0000000000004ebc <+20156>:	movaps 0xdd0(%r8),%xmm15
   0x0000000000004ec4 <+20164>:	mulps  0x150(%rcx),%xmm15
   0x0000000000004ee0 <+20192>:	addps  %xmm10,%xmm15

1513
1514	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1515	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
1516	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004dce <+19918>:	movaps 0xd20(%r8),%xmm10
   0x0000000000004dd6 <+19926>:	mulps  0x160(%rcx),%xmm10
   0x0000000000004df2 <+19954>:	addps  %xmm7,%xmm10
   0x0000000000004e09 <+19977>:	movaps 0xd60(%r8),%xmm10
   0x0000000000004e11 <+19985>:	mulps  0x160(%rcx),%xmm10
   0x0000000000004e45 <+20037>:	movaps 0xda0(%r8),%xmm11
   0x0000000000004e4d <+20045>:	mulps  0x160(%rcx),%xmm11
   0x0000000000004e55 <+20053>:	addps  %xmm8,%xmm10
   0x0000000000004e90 <+20112>:	addps  %xmm9,%xmm11
   0x0000000000004ee4 <+20196>:	movaps 0xde0(%r8),%xmm10
   0x0000000000004eec <+20204>:	mulps  0x160(%rcx),%xmm10
   0x0000000000004f08 <+20232>:	addps  %xmm15,%xmm10

1517
1518	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1519	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
1520	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004d38 <+19768>:	movaps 0xd30(%r8),%xmm14
   0x0000000000004d40 <+19776>:	mulps  0x170(%rcx),%xmm14
   0x0000000000004d82 <+19842>:	movaps 0xdb0(%r8),%xmm12
   0x0000000000004d8a <+19850>:	mulps  0x170(%rcx),%xmm12
   0x0000000000004d96 <+19862>:	movaps 0xd70(%r8),%xmm13
   0x0000000000004d9e <+19870>:	mulps  0x170(%rcx),%xmm13
   0x0000000000004e05 <+19973>:	addps  %xmm10,%xmm14
   0x0000000000004e69 <+20073>:	addps  %xmm10,%xmm13
   0x0000000000004ea4 <+20132>:	addps  %xmm11,%xmm12
   0x0000000000004ea8 <+20136>:	movaps 0xdf0(%r8),%xmm11
   0x0000000000004eb0 <+20144>:	mulps  0x170(%rcx),%xmm11
   0x0000000000004f1c <+20252>:	addps  %xmm10,%xmm11

1521	      }
1522
1523	      /* A(4,3)*B(3,2) = C(4,2). */
1524	      for (i = 0; i < 4; i++)
1525	      {
1526	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1527	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
1528	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004e59 <+20057>:	movaps 0xe40(%r8),%xmm8
   0x0000000000004e61 <+20065>:	mulps  0x240(%rcx),%xmm8
   0x0000000000004e81 <+20097>:	movaps 0xe80(%r8),%xmm7
   0x0000000000004e89 <+20105>:	mulps  0x240(%rcx),%xmm7
   0x0000000000004e94 <+20116>:	movaps 0xe00(%r8),%xmm9
   0x0000000000004e9c <+20124>:	mulps  0x240(%rcx),%xmm9
   0x0000000000004ecc <+20172>:	addps  %xmm14,%xmm9
   0x0000000000004f58 <+20312>:	addps  %xmm13,%xmm8
   0x0000000000004f6c <+20332>:	addps  %xmm12,%xmm7
   0x0000000000004fab <+20395>:	movaps 0xec0(%r8),%xmm13
   0x0000000000004fb3 <+20403>:	mulps  0x240(%rcx),%xmm13
   0x0000000000004ff6 <+20470>:	addps  %xmm11,%xmm13

1529
1530	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1531	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
1532	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004ed0 <+20176>:	movaps 0xe10(%r8),%xmm14
   0x0000000000004ed8 <+20184>:	mulps  0x250(%rcx),%xmm14
   0x0000000000004ef4 <+20212>:	addps  %xmm9,%xmm14
   0x0000000000004f5c <+20316>:	movaps 0xe50(%r8),%xmm13
   0x0000000000004f64 <+20324>:	mulps  0x250(%rcx),%xmm13
   0x0000000000004f70 <+20336>:	movaps 0xe90(%r8),%xmm12
   0x0000000000004f78 <+20344>:	mulps  0x250(%rcx),%xmm12
   0x0000000000004f80 <+20352>:	addps  %xmm8,%xmm13
   0x0000000000004f94 <+20372>:	addps  %xmm7,%xmm12
   0x0000000000004ffa <+20474>:	movaps 0xed0(%r8),%xmm11
   0x0000000000005002 <+20482>:	mulps  0x250(%rcx),%xmm11
   0x000000000000501e <+20510>:	addps  %xmm13,%xmm11

1533
1534	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1535	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
1536	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004ef8 <+20216>:	movaps 0xe20(%r8),%xmm9
   0x0000000000004f00 <+20224>:	mulps  0x260(%rcx),%xmm9
   0x0000000000004f30 <+20272>:	addps  %xmm14,%xmm9
   0x0000000000004f34 <+20276>:	movaps 0xee0(%r8),%xmm14
   0x0000000000004f3c <+20284>:	mulps  0x260(%rcx),%xmm14
   0x0000000000004f84 <+20356>:	movaps 0xe60(%r8),%xmm8
   0x0000000000004f8c <+20364>:	mulps  0x260(%rcx),%xmm8
   0x0000000000004f98 <+20376>:	movaps 0xea0(%r8),%xmm7
   0x0000000000004fa0 <+20384>:	mulps  0x260(%rcx),%xmm7
   0x0000000000004fa7 <+20391>:	addps  %xmm13,%xmm8
   0x0000000000004fcf <+20431>:	addps  %xmm12,%xmm7
   0x0000000000005032 <+20530>:	addps  %xmm11,%xmm14

1537
1538	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1539	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
1540	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004f20 <+20256>:	movaps 0xe30(%r8),%xmm10
   0x0000000000004f28 <+20264>:	mulps  0x270(%rcx),%xmm10
   0x0000000000004f44 <+20292>:	addps  %xmm9,%xmm10
   0x0000000000004f48 <+20296>:	movaps 0xe70(%r8),%xmm9
   0x0000000000004f50 <+20304>:	mulps  0x270(%rcx),%xmm9
   0x0000000000004fbb <+20411>:	addps  %xmm8,%xmm9
   0x0000000000004fbf <+20415>:	movaps 0xeb0(%r8),%xmm8
   0x0000000000004fc7 <+20423>:	mulps  0x270(%rcx),%xmm8
   0x0000000000004fe3 <+20451>:	addps  %xmm7,%xmm8
   0x0000000000004fe7 <+20455>:	movaps 0xef0(%r8),%xmm7
   0x0000000000004fef <+20463>:	mulps  0x270(%rcx),%xmm7
   0x0000000000005046 <+20550>:	addps  %xmm14,%xmm7

1541	      }
1542
1543	      /* A(4,4)*B(4,2) = C(4,2). */
1544	      for (i = 0; i < 4; i++)
1545	      {
1546	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1547	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
1548	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004f0c <+20236>:	movaps 0xf00(%r8),%xmm15
   0x0000000000004f14 <+20244>:	mulps  0x340(%rcx),%xmm15
   0x0000000000004fd3 <+20435>:	movaps 0xf40(%r8),%xmm12
   0x0000000000004fdb <+20443>:	mulps  0x340(%rcx),%xmm12
   0x000000000000500a <+20490>:	addps  %xmm10,%xmm15
   0x000000000000504a <+20554>:	movaps 0xf80(%r8),%xmm14
   0x0000000000005052 <+20562>:	mulps  0x340(%rcx),%xmm14
   0x0000000000005072 <+20594>:	addps  %xmm9,%xmm12
   0x000000000000509e <+20638>:	movaps 0xfc0(%r8),%xmm12
   0x00000000000050a6 <+20646>:	mulps  0x340(%rcx),%xmm12
   0x00000000000050c2 <+20674>:	addps  %xmm8,%xmm14
   0x00000000000050da <+20698>:	addps  %xmm7,%xmm12

1549
1550	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1551	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
1552	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000500e <+20494>:	movaps 0xf10(%r8),%xmm10
   0x0000000000005016 <+20502>:	mulps  0x350(%rcx),%xmm10
   0x000000000000505a <+20570>:	addps  %xmm15,%xmm10
   0x0000000000005076 <+20598>:	movaps 0xf50(%r8),%xmm9
   0x000000000000507e <+20606>:	mulps  0x350(%rcx),%xmm9
   0x000000000000509a <+20634>:	addps  %xmm12,%xmm9
   0x00000000000050c6 <+20678>:	movaps 0xf90(%r8),%xmm8
   0x00000000000050ce <+20686>:	mulps  0x350(%rcx),%xmm8
   0x00000000000050de <+20702>:	movaps 0xfd0(%r8),%xmm7
   0x00000000000050e6 <+20710>:	mulps  0x350(%rcx),%xmm7
   0x00000000000050ed <+20717>:	addps  %xmm14,%xmm8
   0x0000000000005105 <+20741>:	addps  %xmm12,%xmm7

1553
1554	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1555	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
1556	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005022 <+20514>:	movaps 0xf60(%r8),%xmm13
   0x000000000000502a <+20522>:	mulps  0x360(%rcx),%xmm13
   0x0000000000005036 <+20534>:	movaps 0xf20(%r8),%xmm11
   0x000000000000503e <+20542>:	mulps  0x360(%rcx),%xmm11
   0x000000000000505e <+20574>:	addps  %xmm10,%xmm11
   0x000000000000508a <+20618>:	movaps 0xfa0(%r8),%xmm11
   0x0000000000005092 <+20626>:	mulps  0x360(%rcx),%xmm11
   0x00000000000050ae <+20654>:	addps  %xmm9,%xmm13
   0x00000000000050f1 <+20721>:	addps  %xmm8,%xmm11
   0x000000000000511d <+20765>:	movaps 0xfe0(%r8),%xmm11
   0x0000000000005125 <+20773>:	mulps  0x360(%rcx),%xmm11
   0x000000000000512d <+20781>:	addps  %xmm7,%xmm11

1557
1558	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1559	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
1560	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005062 <+20578>:	movaps 0xf30(%r8),%xmm10
   0x000000000000506a <+20586>:	mulps  0x370(%rcx),%xmm10
   0x0000000000005086 <+20614>:	addps  %xmm11,%xmm10
   0x00000000000050b2 <+20658>:	movaps 0xf70(%r8),%xmm9
   0x00000000000050ba <+20666>:	mulps  0x370(%rcx),%xmm9
   0x00000000000050d6 <+20694>:	addps  %xmm13,%xmm9
   0x00000000000050f5 <+20725>:	movaps 0xfb0(%r8),%xmm8
   0x00000000000050fd <+20733>:	mulps  0x370(%rcx),%xmm8
   0x0000000000005109 <+20745>:	movaps 0xff0(%r8),%xmm12
   0x0000000000005111 <+20753>:	mulps  0x370(%rcx),%xmm12
   0x0000000000005119 <+20761>:	addps  %xmm11,%xmm8
   0x0000000000005131 <+20785>:	addps  %xmm11,%xmm12

1561	      }
1562	    }
1563
1564	    /* Store C(4,2) block. */
1565	    for (i = 0; i < 4; i++)
1566	    {
1567	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000005135 <+20789>:	mulps  %xmm0,%xmm10
   0x0000000000005150 <+20816>:	mulps  %xmm0,%xmm9
   0x0000000000005168 <+20840>:	mulps  %xmm0,%xmm8
   0x0000000000005180 <+20864>:	mulps  %xmm0,%xmm12

1568	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_42]), C_row[i]);
   0x0000000000005139 <+20793>:	addps  0x340(%r9),%xmm10
   0x0000000000005154 <+20820>:	addps  0x350(%r9),%xmm9
   0x000000000000516c <+20844>:	addps  0x360(%r9),%xmm8
   0x0000000000005184 <+20868>:	addps  0x370(%r9),%xmm12

1569	      _mm_store_ps(&C[i*4+C_OFFSET_42], C_row[i]);
   0x0000000000005148 <+20808>:	movaps %xmm10,0x340(%r9)
   0x000000000000515c <+20828>:	movaps %xmm9,0x350(%r9)
   0x0000000000005174 <+20852>:	movaps %xmm8,0x360(%r9)
   0x000000000000518c <+20876>:	movaps %xmm12,0x370(%r9)

1570	    }
1571
1572	    /* Reset C(4,3) matrix accumulators */
1573	    C_row[0] = _mm_setzero_ps();
   0x0000000000005164 <+20836>:	xorps  %xmm9,%xmm9

1574	    C_row[1] = _mm_setzero_ps();
   0x000000000000517c <+20860>:	xorps  %xmm8,%xmm8

1575	    C_row[2] = _mm_setzero_ps();
   0x0000000000005141 <+20801>:	xorps  %xmm7,%xmm7

1576	    C_row[3] = _mm_setzero_ps();
   0x0000000000005144 <+20804>:	xorps  %xmm11,%xmm11

1577
1578	    if (norm_product[12][2] &&
   0x0000000000005194 <+20884>:	comiss %xmm1,%xmm6
   0x0000000000005197 <+20887>:	jb     0x56a4 <stream_kernel+22180>

1579	        norm_product[13][6] &&
   0x000000000000519d <+20893>:	movss  0x58(%rsp),%xmm6
   0x00000000000051a3 <+20899>:	comiss %xmm1,%xmm6
   0x00000000000051a6 <+20902>:	jb     0x56a4 <stream_kernel+22180>

1580	        norm_product[14][10] &&
   0x00000000000051ac <+20908>:	movss  0x20(%rsp),%xmm6
   0x00000000000051b2 <+20914>:	comiss %xmm1,%xmm6
   0x00000000000051b5 <+20917>:	jb     0x56a4 <stream_kernel+22180>

1581	        norm_product[15][14])
   0x00000000000051bb <+20923>:	movss  0x10(%rsp),%xmm6
   0x00000000000051c1 <+20929>:	comiss %xmm1,%xmm6
   0x00000000000051c4 <+20932>:	jb     0x56a4 <stream_kernel+22180>

1582	    {
1583	      /* A(4,1)*B(1,3) = C(4,3). */
1584	      for (i = 0; i < 4; i++)
1585	      {
1586	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1587	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
1588	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000051ca <+20938>:	movaps 0xc00(%r8),%xmm7
   0x00000000000051ea <+20970>:	movaps 0xc40(%r8),%xmm11
   0x000000000000520a <+21002>:	movaps 0xc80(%r8),%xmm14
   0x000000000000521a <+21018>:	mulps  0x80(%rcx),%xmm7
   0x0000000000005238 <+21048>:	mulps  0x80(%rcx),%xmm11
   0x0000000000005258 <+21080>:	mulps  0x80(%rcx),%xmm14
   0x0000000000005291 <+21137>:	movaps 0xcc0(%r8),%xmm10
   0x0000000000005299 <+21145>:	mulps  0x80(%rcx),%xmm10

1589
1590	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1591	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
1592	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000051d2 <+20946>:	movaps 0xc10(%r8),%xmm6
   0x00000000000051f2 <+20978>:	movaps 0xc50(%r8),%xmm12
   0x0000000000005221 <+21025>:	mulps  0x90(%rcx),%xmm6
   0x0000000000005240 <+21056>:	mulps  0x90(%rcx),%xmm12
   0x0000000000005268 <+21096>:	addps  %xmm7,%xmm6
   0x000000000000526b <+21099>:	movaps 0xc90(%r8),%xmm7
   0x0000000000005273 <+21107>:	mulps  0x90(%rcx),%xmm7
   0x00000000000052a1 <+21153>:	addps  %xmm11,%xmm12
   0x00000000000052a5 <+21157>:	movaps 0xcd0(%r8),%xmm11
   0x00000000000052ad <+21165>:	mulps  0x90(%rcx),%xmm11
   0x00000000000052dd <+21213>:	addps  %xmm14,%xmm7
   0x0000000000005316 <+21270>:	addps  %xmm10,%xmm11

1593
1594	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1595	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
1596	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000051da <+20954>:	movaps 0xc20(%r8),%xmm10
   0x00000000000051fa <+20986>:	movaps 0xc60(%r8),%xmm13
   0x0000000000005228 <+21032>:	mulps  0xa0(%rcx),%xmm10
   0x0000000000005248 <+21064>:	mulps  0xa0(%rcx),%xmm13
   0x000000000000527a <+21114>:	addps  %xmm6,%xmm10
   0x000000000000527e <+21118>:	movaps 0xca0(%r8),%xmm6
   0x0000000000005286 <+21126>:	mulps  0xa0(%rcx),%xmm6
   0x00000000000052b5 <+21173>:	addps  %xmm12,%xmm13
   0x00000000000052b9 <+21177>:	movaps 0xce0(%r8),%xmm12
   0x00000000000052c1 <+21185>:	mulps  0xa0(%rcx),%xmm12
   0x00000000000052f1 <+21233>:	addps  %xmm7,%xmm6
   0x000000000000532a <+21290>:	addps  %xmm11,%xmm12

1597
1598	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1599	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
1600	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000051e2 <+20962>:	movaps 0xc30(%r8),%xmm8
   0x0000000000005202 <+20994>:	movaps 0xc70(%r8),%xmm9
   0x0000000000005212 <+21010>:	movaps 0xcb0(%r8),%xmm15
   0x0000000000005230 <+21040>:	mulps  0xb0(%rcx),%xmm8
   0x0000000000005250 <+21072>:	mulps  0xb0(%rcx),%xmm9
   0x0000000000005260 <+21088>:	mulps  0xb0(%rcx),%xmm15
   0x000000000000528d <+21133>:	addps  %xmm10,%xmm8
   0x00000000000052c9 <+21193>:	addps  %xmm13,%xmm9
   0x00000000000052e1 <+21217>:	movaps 0xcf0(%r8),%xmm14
   0x00000000000052e9 <+21225>:	mulps  0xb0(%rcx),%xmm14
   0x0000000000005303 <+21251>:	addps  %xmm6,%xmm15
   0x000000000000533e <+21310>:	addps  %xmm12,%xmm14

1601	      }
1602
1603	      /* A(4,2)*B(2,3) = C(4,3). */
1604	      for (i = 0; i < 4; i++)
1605	      {
1606	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1607	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
1608	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000052f4 <+21236>:	movaps 0xd00(%r8),%xmm7
   0x00000000000052fc <+21244>:	mulps  0x180(%rcx),%xmm7
   0x0000000000005307 <+21255>:	movaps 0xd40(%r8),%xmm6
   0x000000000000530f <+21263>:	mulps  0x180(%rcx),%xmm6
   0x0000000000005352 <+21330>:	addps  %xmm8,%xmm7
   0x0000000000005366 <+21350>:	addps  %xmm9,%xmm6
   0x00000000000053b8 <+21432>:	movaps 0xd80(%r8),%xmm7
   0x00000000000053c0 <+21440>:	mulps  0x180(%rcx),%xmm7
   0x00000000000053df <+21471>:	movaps 0xdc0(%r8),%xmm8
   0x00000000000053e7 <+21479>:	mulps  0x180(%rcx),%xmm8
   0x00000000000053ef <+21487>:	addps  %xmm15,%xmm7
   0x0000000000005403 <+21507>:	addps  %xmm14,%xmm8

1609
1610	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1611	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
1612	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005356 <+21334>:	movaps 0xd10(%r8),%xmm8
   0x000000000000535e <+21342>:	mulps  0x190(%rcx),%xmm8
   0x000000000000536a <+21354>:	movaps 0xd50(%r8),%xmm9
   0x0000000000005372 <+21362>:	mulps  0x190(%rcx),%xmm9
   0x000000000000537a <+21370>:	addps  %xmm7,%xmm8
   0x000000000000538d <+21389>:	addps  %xmm6,%xmm9
   0x00000000000053f3 <+21491>:	movaps 0xd90(%r8),%xmm15
   0x00000000000053fb <+21499>:	mulps  0x190(%rcx),%xmm15
   0x0000000000005417 <+21527>:	addps  %xmm7,%xmm15
   0x000000000000541b <+21531>:	movaps 0xdd0(%r8),%xmm7
   0x0000000000005423 <+21539>:	mulps  0x190(%rcx),%xmm7
   0x0000000000005452 <+21586>:	addps  %xmm8,%xmm7

1613
1614	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1615	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
1616	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000537e <+21374>:	movaps 0xd20(%r8),%xmm7
   0x0000000000005386 <+21382>:	mulps  0x1a0(%rcx),%xmm7
   0x0000000000005391 <+21393>:	movaps 0xde0(%r8),%xmm6
   0x0000000000005399 <+21401>:	mulps  0x1a0(%rcx),%xmm6
   0x00000000000053a0 <+21408>:	addps  %xmm8,%xmm7
   0x00000000000053a4 <+21412>:	movaps 0xd60(%r8),%xmm8
   0x00000000000053ac <+21420>:	mulps  0x1a0(%rcx),%xmm8
   0x00000000000053c7 <+21447>:	addps  %xmm9,%xmm8
   0x00000000000053cb <+21451>:	movaps 0xda0(%r8),%xmm9
   0x00000000000053d3 <+21459>:	mulps  0x1a0(%rcx),%xmm9
   0x000000000000542a <+21546>:	addps  %xmm15,%xmm9
   0x0000000000005466 <+21606>:	addps  %xmm7,%xmm6

1617
1618	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1619	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
1620	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000052cd <+21197>:	movaps 0xd30(%r8),%xmm13
   0x00000000000052d5 <+21205>:	mulps  0x1b0(%rcx),%xmm13
   0x000000000000531a <+21274>:	movaps 0xdf0(%r8),%xmm10
   0x0000000000005322 <+21282>:	mulps  0x1b0(%rcx),%xmm10
   0x000000000000532e <+21294>:	movaps 0xdb0(%r8),%xmm11
   0x0000000000005336 <+21302>:	mulps  0x1b0(%rcx),%xmm11
   0x0000000000005342 <+21314>:	movaps 0xd70(%r8),%xmm12
   0x000000000000534a <+21322>:	mulps  0x1b0(%rcx),%xmm12
   0x00000000000053b4 <+21428>:	addps  %xmm7,%xmm13
   0x00000000000053db <+21467>:	addps  %xmm8,%xmm12
   0x000000000000543e <+21566>:	addps  %xmm9,%xmm11
   0x0000000000005478 <+21624>:	addps  %xmm6,%xmm10

1621	      }
1622
1623	      /* A(4,3)*B(3,3) = C(4,3). */
1624	      for (i = 0; i < 4; i++)
1625	      {
1626	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1627	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
1628	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005407 <+21511>:	movaps 0xe00(%r8),%xmm14
   0x000000000000540f <+21519>:	mulps  0x280(%rcx),%xmm14
   0x0000000000005469 <+21609>:	movaps 0xe40(%r8),%xmm7
   0x0000000000005471 <+21617>:	mulps  0x280(%rcx),%xmm7
   0x000000000000548b <+21643>:	addps  %xmm13,%xmm14
   0x000000000000549f <+21663>:	addps  %xmm12,%xmm7
   0x00000000000054b7 <+21687>:	movaps 0xe80(%r8),%xmm14
   0x00000000000054bf <+21695>:	mulps  0x280(%rcx),%xmm14
   0x0000000000005505 <+21765>:	movaps 0xec0(%r8),%xmm6
   0x000000000000550d <+21773>:	mulps  0x280(%rcx),%xmm6
   0x0000000000005528 <+21800>:	addps  %xmm11,%xmm14
   0x000000000000553c <+21820>:	addps  %xmm10,%xmm6

1629
1630	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1631	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
1632	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000547c <+21628>:	movaps 0xe50(%r8),%xmm6
   0x0000000000005484 <+21636>:	mulps  0x290(%rcx),%xmm6
   0x000000000000548f <+21647>:	movaps 0xe10(%r8),%xmm13
   0x0000000000005497 <+21655>:	mulps  0x290(%rcx),%xmm13
   0x00000000000054b3 <+21683>:	addps  %xmm14,%xmm13
   0x00000000000054ef <+21743>:	addps  %xmm7,%xmm6
   0x000000000000552c <+21804>:	movaps 0xe90(%r8),%xmm11
   0x0000000000005534 <+21812>:	mulps  0x290(%rcx),%xmm11
   0x0000000000005540 <+21824>:	movaps 0xed0(%r8),%xmm10
   0x0000000000005548 <+21832>:	mulps  0x290(%rcx),%xmm10
   0x0000000000005550 <+21840>:	addps  %xmm14,%xmm11
   0x000000000000558c <+21900>:	addps  %xmm6,%xmm10

1633
1634	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1635	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
1636	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000542e <+21550>:	movaps 0xea0(%r8),%xmm15
   0x0000000000005436 <+21558>:	mulps  0x2a0(%rcx),%xmm15
   0x0000000000005456 <+21590>:	movaps 0xe20(%r8),%xmm8
   0x000000000000545e <+21598>:	mulps  0x2a0(%rcx),%xmm8
   0x00000000000054a3 <+21667>:	movaps 0xe60(%r8),%xmm12
   0x00000000000054ab <+21675>:	mulps  0x2a0(%rcx),%xmm12
   0x00000000000054c7 <+21703>:	addps  %xmm13,%xmm8
   0x0000000000005501 <+21761>:	addps  %xmm6,%xmm12
   0x0000000000005564 <+21860>:	addps  %xmm11,%xmm15
   0x0000000000005568 <+21864>:	movaps 0xee0(%r8),%xmm11
   0x0000000000005570 <+21872>:	mulps  0x2a0(%rcx),%xmm11
   0x000000000000559f <+21919>:	addps  %xmm10,%xmm11

1637
1638	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1639	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
1640	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005442 <+21570>:	movaps 0xe30(%r8),%xmm9
   0x000000000000544a <+21578>:	mulps  0x2b0(%rcx),%xmm9
   0x00000000000054db <+21723>:	addps  %xmm8,%xmm9
   0x00000000000054df <+21727>:	movaps 0xe70(%r8),%xmm8
   0x00000000000054e7 <+21735>:	mulps  0x2b0(%rcx),%xmm8
   0x00000000000054f2 <+21746>:	movaps 0xeb0(%r8),%xmm7
   0x00000000000054fa <+21754>:	mulps  0x2b0(%rcx),%xmm7
   0x0000000000005514 <+21780>:	addps  %xmm12,%xmm8
   0x0000000000005578 <+21880>:	addps  %xmm15,%xmm7
   0x0000000000005590 <+21904>:	movaps 0xef0(%r8),%xmm6
   0x0000000000005598 <+21912>:	mulps  0x2b0(%rcx),%xmm6
   0x00000000000055b3 <+21939>:	addps  %xmm11,%xmm6

1641	      }
1642
1643	      /* A(4,4)*B(4,3) = C(4,3). */
1644	      for (i = 0; i < 4; i++)
1645	      {
1646	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1647	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
1648	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005518 <+21784>:	movaps 0xf00(%r8),%xmm12
   0x0000000000005520 <+21792>:	mulps  0x380(%rcx),%xmm12
   0x0000000000005554 <+21844>:	movaps 0xf40(%r8),%xmm14
   0x000000000000555c <+21852>:	mulps  0x380(%rcx),%xmm14
   0x00000000000055a3 <+21923>:	movaps 0xf80(%r8),%xmm10
   0x00000000000055ab <+21931>:	mulps  0x380(%rcx),%xmm10
   0x00000000000055c7 <+21959>:	addps  %xmm9,%xmm12
   0x00000000000055db <+21979>:	addps  %xmm8,%xmm14
   0x00000000000055f3 <+22003>:	movaps 0xfc0(%r8),%xmm12
   0x00000000000055fb <+22011>:	mulps  0x380(%rcx),%xmm12
   0x0000000000005633 <+22067>:	addps  %xmm7,%xmm10
   0x000000000000564a <+22090>:	addps  %xmm6,%xmm12

1649
1650	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1651	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
1652	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000055cb <+21963>:	movaps 0xf10(%r8),%xmm9
   0x00000000000055d3 <+21971>:	mulps  0x390(%rcx),%xmm9
   0x00000000000055df <+21983>:	movaps 0xf50(%r8),%xmm8
   0x00000000000055e7 <+21991>:	mulps  0x390(%rcx),%xmm8
   0x00000000000055ef <+21999>:	addps  %xmm12,%xmm9
   0x0000000000005617 <+22039>:	addps  %xmm14,%xmm8
   0x0000000000005637 <+22071>:	movaps 0xf90(%r8),%xmm7
   0x000000000000563f <+22079>:	mulps  0x390(%rcx),%xmm7
   0x000000000000564e <+22094>:	movaps 0xfd0(%r8),%xmm6
   0x0000000000005656 <+22102>:	mulps  0x390(%rcx),%xmm6
   0x000000000000565d <+22109>:	addps  %xmm10,%xmm7
   0x0000000000005684 <+22148>:	addps  %xmm12,%xmm6

1653
1654	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1655	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
1656	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000054cb <+21707>:	movaps 0xf20(%r8),%xmm13
   0x00000000000054d3 <+21715>:	mulps  0x3a0(%rcx),%xmm13
   0x000000000000557c <+21884>:	movaps 0xf60(%r8),%xmm15
   0x0000000000005584 <+21892>:	mulps  0x3a0(%rcx),%xmm15
   0x00000000000055b7 <+21943>:	movaps 0xfa0(%r8),%xmm11
   0x00000000000055bf <+21951>:	mulps  0x3a0(%rcx),%xmm11
   0x0000000000005603 <+22019>:	addps  %xmm9,%xmm13
   0x000000000000561b <+22043>:	addps  %xmm8,%xmm15
   0x0000000000005661 <+22113>:	movaps 0xfe0(%r8),%xmm10
   0x0000000000005669 <+22121>:	mulps  0x3a0(%rcx),%xmm10
   0x0000000000005671 <+22129>:	addps  %xmm7,%xmm11
   0x0000000000005688 <+22152>:	addps  %xmm6,%xmm10

1657
1658	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1659	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
1660	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005607 <+22023>:	movaps 0xf30(%r8),%xmm9
   0x000000000000560f <+22031>:	mulps  0x3b0(%rcx),%xmm9
   0x000000000000561f <+22047>:	movaps 0xf70(%r8),%xmm8
   0x0000000000005627 <+22055>:	mulps  0x3b0(%rcx),%xmm8
   0x000000000000562f <+22063>:	addps  %xmm13,%xmm9
   0x0000000000005646 <+22086>:	addps  %xmm15,%xmm8
   0x0000000000005675 <+22133>:	movaps 0xfb0(%r8),%xmm7
   0x000000000000567d <+22141>:	mulps  0x3b0(%rcx),%xmm7
   0x000000000000568c <+22156>:	addps  %xmm11,%xmm7
   0x0000000000005690 <+22160>:	movaps 0xff0(%r8),%xmm11
   0x0000000000005698 <+22168>:	mulps  0x3b0(%rcx),%xmm11
   0x00000000000056a0 <+22176>:	addps  %xmm10,%xmm11

1661	      }
1662	    }
1663
1664	    /* Store C(4,3) block. */
1665	    for (i = 0; i < 4; i++)
1666	    {
1667	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000056a4 <+22180>:	mulps  %xmm0,%xmm9
   0x00000000000056c3 <+22211>:	mulps  %xmm0,%xmm8
   0x00000000000056db <+22235>:	mulps  %xmm0,%xmm7
   0x00000000000056ee <+22254>:	mulps  %xmm0,%xmm11

1668	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_43]), C_row[i]);
   0x00000000000056a8 <+22184>:	addps  0x380(%r9),%xmm9
   0x00000000000056c7 <+22215>:	addps  0x390(%r9),%xmm8
   0x00000000000056de <+22238>:	addps  0x3a0(%r9),%xmm7
   0x00000000000056f2 <+22258>:	addps  0x3b0(%r9),%xmm11

1669	      _mm_store_ps(&C[i*4+C_OFFSET_43], C_row[i]);
   0x00000000000056b7 <+22199>:	movaps %xmm9,0x380(%r9)
   0x00000000000056cf <+22223>:	movaps %xmm8,0x390(%r9)
   0x00000000000056e6 <+22246>:	movaps %xmm7,0x3a0(%r9)
   0x00000000000056fa <+22266>:	movaps %xmm11,0x3b0(%r9)

1670	    }
1671
1672	    /* Reset C(4,4) matrix accumulators */
1673	    C_row[0] = _mm_setzero_ps();
   0x00000000000056b0 <+22192>:	xorps  %xmm12,%xmm12

1674	    C_row[1] = _mm_setzero_ps();
   0x00000000000056d7 <+22231>:	xorps  %xmm8,%xmm8

1675	    C_row[2] = _mm_setzero_ps();
   0x00000000000056b4 <+22196>:	xorps  %xmm6,%xmm6

1676	    C_row[3] = _mm_setzero_ps();
   0x00000000000056bf <+22207>:	xorps  %xmm9,%xmm9

1677
1678	    if (norm_product[12][3] &&
   0x0000000000005702 <+22274>:	comiss %xmm1,%xmm5
   0x0000000000005705 <+22277>:	jb     0x5b7d <stream_kernel+23421>

1679	        norm_product[13][7] &&
   0x000000000000570b <+22283>:	comiss %xmm1,%xmm4
   0x000000000000570e <+22286>:	jb     0x5b7d <stream_kernel+23421>

1680	        norm_product[14][11] &&
   0x0000000000005714 <+22292>:	comiss %xmm1,%xmm3
   0x0000000000005717 <+22295>:	jb     0x5b7d <stream_kernel+23421>

1681	        norm_product[15][15])
   0x000000000000571d <+22301>:	comiss %xmm1,%xmm2
   0x0000000000005720 <+22304>:	jb     0x5b7d <stream_kernel+23421>

1682	    {
1683	      /* A(4,1)*B(1,4) = C(4,4). */
1684	      for (i = 0; i < 4; i++)
1685	      {
1686	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1687	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
1688	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005726 <+22310>:	movaps 0xc0(%rcx),%xmm15
   0x0000000000005745 <+22341>:	movaps 0xc00(%r8),%xmm9
   0x0000000000005765 <+22373>:	movaps 0xc40(%r8),%xmm8
   0x0000000000005785 <+22405>:	movaps 0xc80(%r8),%xmm13
   0x0000000000005795 <+22421>:	mulps  %xmm15,%xmm9
   0x00000000000057b1 <+22449>:	movaps 0xcc0(%r8),%xmm12
   0x00000000000057c6 <+22470>:	mulps  %xmm15,%xmm8
   0x00000000000057f9 <+22521>:	mulps  %xmm15,%xmm13
   0x000000000000582c <+22572>:	mulps  %xmm15,%xmm12

1689
1690	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1691	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
1692	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000572e <+22318>:	movaps 0xd0(%rcx),%xmm14
   0x000000000000574d <+22349>:	movaps 0xc10(%r8),%xmm12
   0x000000000000576d <+22381>:	movaps 0xc50(%r8),%xmm6
   0x000000000000578d <+22413>:	movaps 0xc90(%r8),%xmm4
   0x0000000000005799 <+22425>:	mulps  %xmm14,%xmm12
   0x000000000000579d <+22429>:	addps  %xmm9,%xmm12
   0x00000000000057ca <+22474>:	mulps  %xmm14,%xmm6
   0x00000000000057ce <+22478>:	addps  %xmm8,%xmm6
   0x00000000000057fd <+22525>:	mulps  %xmm14,%xmm4
   0x0000000000005801 <+22529>:	addps  %xmm13,%xmm4
   0x0000000000005830 <+22576>:	movaps 0xcd0(%r8),%xmm15
   0x0000000000005838 <+22584>:	mulps  %xmm14,%xmm15
   0x0000000000005844 <+22596>:	addps  %xmm12,%xmm15

1693
1694	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1695	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
1696	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005736 <+22326>:	movaps 0xe0(%rcx),%xmm11
   0x0000000000005755 <+22357>:	movaps 0xc20(%r8),%xmm7
   0x0000000000005775 <+22389>:	movaps 0xc60(%r8),%xmm10
   0x00000000000057a1 <+22433>:	movaps 0xca0(%r8),%xmm9
   0x00000000000057a9 <+22441>:	mulps  %xmm11,%xmm7
   0x00000000000057ad <+22445>:	addps  %xmm12,%xmm7
   0x00000000000057da <+22490>:	mulps  %xmm11,%xmm10
   0x00000000000057de <+22494>:	addps  %xmm6,%xmm10
   0x000000000000580d <+22541>:	mulps  %xmm11,%xmm9
   0x0000000000005811 <+22545>:	addps  %xmm4,%xmm9
   0x000000000000583c <+22588>:	movaps 0xce0(%r8),%xmm14
   0x0000000000005848 <+22600>:	mulps  %xmm11,%xmm14
   0x0000000000005854 <+22612>:	addps  %xmm15,%xmm14

1697
1698	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1699	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
1700	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000573e <+22334>:	movaps 0xf0(%rcx),%xmm5
   0x000000000000575d <+22365>:	movaps 0xc30(%r8),%xmm3
   0x000000000000577d <+22397>:	movaps 0xc70(%r8),%xmm2
   0x00000000000057b9 <+22457>:	mulps  %xmm5,%xmm3
   0x00000000000057bc <+22460>:	addps  %xmm7,%xmm3
   0x00000000000057ea <+22506>:	mulps  %xmm5,%xmm2
   0x00000000000057ed <+22509>:	addps  %xmm10,%xmm2
   0x0000000000005805 <+22533>:	movaps 0xcf0(%r8),%xmm13
   0x0000000000005815 <+22549>:	movaps 0xcb0(%r8),%xmm4
   0x000000000000581d <+22557>:	mulps  %xmm5,%xmm4
   0x0000000000005820 <+22560>:	addps  %xmm9,%xmm4
   0x0000000000005858 <+22616>:	mulps  %xmm5,%xmm13
   0x0000000000005864 <+22628>:	addps  %xmm14,%xmm13

1701	      }
1702
1703	      /* A(4,2)*B(2,4) = C(4,4). */
1704	      for (i = 0; i < 4; i++)
1705	      {
1706	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1707	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
1708	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000057f1 <+22513>:	movaps 0x1c0(%rcx),%xmm10
   0x0000000000005824 <+22564>:	movaps 0xd40(%r8),%xmm9
   0x000000000000584c <+22604>:	movaps 0xd00(%r8),%xmm11
   0x0000000000005868 <+22632>:	mulps  %xmm10,%xmm11
   0x0000000000005874 <+22644>:	addps  %xmm3,%xmm11
   0x000000000000588c <+22668>:	movaps 0xd80(%r8),%xmm15
   0x00000000000058aa <+22698>:	mulps  %xmm10,%xmm9
   0x00000000000058c5 <+22725>:	addps  %xmm2,%xmm9
   0x00000000000058e8 <+22760>:	movaps 0xdc0(%r8),%xmm11
   0x00000000000058f0 <+22768>:	mulps  %xmm10,%xmm15
   0x0000000000005904 <+22788>:	addps  %xmm4,%xmm15
   0x0000000000005923 <+22819>:	mulps  %xmm10,%xmm11
   0x0000000000005967 <+22887>:	addps  %xmm13,%xmm11

1709
1710	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1711	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
1712	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000057d2 <+22482>:	movaps 0x1d0(%rcx),%xmm8
   0x000000000000585c <+22620>:	movaps 0xd10(%r8),%xmm5
   0x0000000000005878 <+22648>:	mulps  %xmm8,%xmm5
   0x0000000000005898 <+22680>:	addps  %xmm11,%xmm5
   0x000000000000589c <+22684>:	movaps 0xd50(%r8),%xmm11
   0x00000000000058b9 <+22713>:	mulps  %xmm8,%xmm11
   0x00000000000058d5 <+22741>:	addps  %xmm9,%xmm11
   0x000000000000590b <+22795>:	movaps 0xd90(%r8),%xmm4
   0x0000000000005913 <+22803>:	mulps  %xmm8,%xmm4
   0x0000000000005917 <+22807>:	addps  %xmm15,%xmm4
   0x0000000000005927 <+22823>:	movaps 0xdd0(%r8),%xmm10
   0x0000000000005933 <+22835>:	mulps  %xmm8,%xmm10
   0x0000000000005973 <+22899>:	addps  %xmm11,%xmm10

1713
1714	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1715	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
1716	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000586c <+22636>:	movaps 0x1e0(%rcx),%xmm12
   0x000000000000587c <+22652>:	movaps 0xd20(%r8),%xmm3
   0x0000000000005884 <+22660>:	movaps 0xd60(%r8),%xmm14
   0x0000000000005894 <+22676>:	mulps  %xmm12,%xmm3
   0x00000000000058a7 <+22695>:	addps  %xmm5,%xmm3
   0x00000000000058d1 <+22737>:	mulps  %xmm12,%xmm14
   0x00000000000058d9 <+22745>:	movaps 0xda0(%r8),%xmm9
   0x00000000000058e4 <+22756>:	addps  %xmm11,%xmm14
   0x00000000000058f8 <+22776>:	mulps  %xmm12,%xmm9
   0x000000000000592f <+22831>:	addps  %xmm4,%xmm9
   0x0000000000005937 <+22839>:	movaps 0xde0(%r8),%xmm8
   0x000000000000594b <+22859>:	mulps  %xmm12,%xmm8
   0x000000000000597f <+22911>:	addps  %xmm10,%xmm8

1717
1718	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1719	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
1720	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000057bf <+22463>:	movaps 0x1f0(%rcx),%xmm7
   0x00000000000057e2 <+22498>:	movaps 0xd30(%r8),%xmm6
   0x00000000000058a4 <+22692>:	mulps  %xmm7,%xmm6
   0x00000000000058b6 <+22710>:	addps  %xmm3,%xmm6
   0x00000000000058bd <+22717>:	movaps 0xdb0(%r8),%xmm3
   0x00000000000058c9 <+22729>:	movaps 0xd70(%r8),%xmm2
   0x00000000000058e1 <+22753>:	mulps  %xmm7,%xmm2
   0x00000000000058f4 <+22772>:	addps  %xmm14,%xmm2
   0x0000000000005908 <+22792>:	mulps  %xmm7,%xmm3
   0x000000000000593f <+22847>:	addps  %xmm9,%xmm3
   0x0000000000005943 <+22851>:	movaps 0xdf0(%r8),%xmm9
   0x0000000000005957 <+22871>:	mulps  %xmm7,%xmm9
   0x0000000000005997 <+22935>:	addps  %xmm8,%xmm9

1721	      }
1722
1723	      /* A(4,3)*B(3,4) = C(4,4). */
1724	      for (i = 0; i < 4; i++)
1725	      {
1726	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1727	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
1728	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000591b <+22811>:	movaps 0x2c0(%rcx),%xmm15
   0x000000000000594f <+22863>:	movaps 0xe00(%r8),%xmm12
   0x0000000000005963 <+22883>:	mulps  %xmm15,%xmm12
   0x0000000000005983 <+22915>:	movaps 0xe80(%r8),%xmm10
   0x000000000000598b <+22923>:	movaps 0xec0(%r8),%xmm11
   0x0000000000005993 <+22931>:	mulps  %xmm15,%xmm10
   0x000000000000599b <+22939>:	mulps  %xmm15,%xmm11
   0x00000000000059a7 <+22951>:	addps  %xmm6,%xmm12
   0x00000000000059b5 <+22965>:	addps  %xmm3,%xmm10
   0x00000000000059d5 <+22997>:	addps  %xmm9,%xmm11
   0x00000000000059ec <+23020>:	movaps 0xe40(%r8),%xmm7
   0x00000000000059f8 <+23032>:	mulps  %xmm15,%xmm7
   0x0000000000005a14 <+23060>:	addps  %xmm2,%xmm7

1729
1730	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1731	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
1732	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000058ae <+22702>:	movaps 0xe90(%r8),%xmm5
   0x000000000000595b <+22875>:	movaps 0xe10(%r8),%xmm7
   0x00000000000059ab <+22955>:	movaps 0x2d0(%rcx),%xmm6
   0x00000000000059b2 <+22962>:	mulps  %xmm6,%xmm7
   0x00000000000059c5 <+22981>:	addps  %xmm12,%xmm7
   0x00000000000059dd <+23005>:	movaps 0xed0(%r8),%xmm9
   0x00000000000059e9 <+23017>:	mulps  %xmm6,%xmm5
   0x0000000000005a04 <+23044>:	addps  %xmm10,%xmm5
   0x0000000000005a17 <+23063>:	mulps  %xmm6,%xmm9
   0x0000000000005a1b <+23067>:	movaps 0xe50(%r8),%xmm2
   0x0000000000005a27 <+23079>:	mulps  %xmm6,%xmm2
   0x0000000000005a3e <+23102>:	addps  %xmm7,%xmm2
   0x0000000000005a8f <+23183>:	addps  %xmm11,%xmm9

1733
1734	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1735	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
1736	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000596b <+22891>:	movaps 0xe20(%r8),%xmm13
   0x00000000000059cd <+22989>:	movaps 0x2e0(%rcx),%xmm12
   0x00000000000059d9 <+23001>:	mulps  %xmm12,%xmm13
   0x00000000000059e5 <+23013>:	addps  %xmm7,%xmm13
   0x00000000000059fc <+23036>:	movaps 0xea0(%r8),%xmm13
   0x0000000000005a08 <+23048>:	mulps  %xmm12,%xmm13
   0x0000000000005a23 <+23075>:	addps  %xmm5,%xmm13
   0x0000000000005a32 <+23090>:	movaps 0xee0(%r8),%xmm6
   0x0000000000005a41 <+23105>:	mulps  %xmm12,%xmm6
   0x0000000000005a45 <+23109>:	movaps 0xe60(%r8),%xmm7
   0x0000000000005a51 <+23121>:	mulps  %xmm12,%xmm7
   0x0000000000005a75 <+23157>:	addps  %xmm2,%xmm7
   0x0000000000005a9f <+23199>:	addps  %xmm9,%xmm6

1737
1738	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1739	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
1740	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000058fc <+22780>:	movaps 0xe30(%r8),%xmm14
   0x0000000000005977 <+22903>:	movaps 0xe70(%r8),%xmm4
   0x000000000000599f <+22943>:	movaps 0x2f0(%rcx),%xmm8
   0x00000000000059b9 <+22969>:	mulps  %xmm8,%xmm14
   0x00000000000059c9 <+22985>:	mulps  %xmm8,%xmm4
   0x00000000000059f4 <+23028>:	addps  %xmm13,%xmm14
   0x0000000000005a2a <+23082>:	movaps 0xeb0(%r8),%xmm5
   0x0000000000005a3a <+23098>:	mulps  %xmm8,%xmm5
   0x0000000000005a4d <+23117>:	addps  %xmm13,%xmm5
   0x0000000000005a55 <+23125>:	movaps 0xef0(%r8),%xmm13
   0x0000000000005a65 <+23141>:	mulps  %xmm8,%xmm13
   0x0000000000005a80 <+23168>:	addps  %xmm7,%xmm4
   0x0000000000005aaf <+23215>:	addps  %xmm6,%xmm13

1741	      }
1742
1743	      /* A(4,4)*B(4,4) = C(4,4). */
1744	      for (i = 0; i < 4; i++)
1745	      {
1746	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1747	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
1748	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005a78 <+23160>:	movaps 0xf40(%r8),%xmm2
   0x0000000000005a93 <+23187>:	movaps 0x3c0(%rcx),%xmm11
   0x0000000000005a9b <+23195>:	mulps  %xmm11,%xmm2
   0x0000000000005aa3 <+23203>:	movaps 0xf00(%r8),%xmm9
   0x0000000000005aab <+23211>:	mulps  %xmm11,%xmm9
   0x0000000000005ac3 <+23235>:	addps  %xmm14,%xmm9
   0x0000000000005ac7 <+23239>:	movaps 0xf80(%r8),%xmm14
   0x0000000000005acf <+23247>:	mulps  %xmm11,%xmm14
   0x0000000000005ae3 <+23267>:	addps  %xmm4,%xmm2
   0x0000000000005b08 <+23304>:	movaps 0xfc0(%r8),%xmm2
   0x0000000000005b22 <+23330>:	mulps  %xmm11,%xmm2
   0x0000000000005b2e <+23342>:	addps  %xmm5,%xmm14
   0x0000000000005b3e <+23358>:	addps  %xmm13,%xmm2

1749
1750	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1751	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
1752	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000059bd <+22973>:	movaps 0xf50(%r8),%xmm3
   0x0000000000005a0c <+23052>:	movaps 0x3d0(%rcx),%xmm10
   0x0000000000005a71 <+23153>:	mulps  %xmm10,%xmm3
   0x0000000000005a83 <+23171>:	movaps 0xf10(%r8),%xmm7
   0x0000000000005a8b <+23179>:	mulps  %xmm10,%xmm7
   0x0000000000005ad3 <+23251>:	addps  %xmm9,%xmm7
   0x0000000000005b01 <+23297>:	addps  %xmm2,%xmm3
   0x0000000000005b13 <+23315>:	movaps 0xfd0(%r8),%xmm3
   0x0000000000005b32 <+23346>:	movaps 0xf90(%r8),%xmm5
   0x0000000000005b42 <+23362>:	mulps  %xmm10,%xmm5
   0x0000000000005b4a <+23370>:	mulps  %xmm10,%xmm3
   0x0000000000005b56 <+23382>:	addps  %xmm14,%xmm5
   0x0000000000005b72 <+23410>:	addps  %xmm2,%xmm3

1753
1754	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1755	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
1756	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005a69 <+23145>:	movaps 0xf20(%r8),%xmm8
   0x0000000000005ab3 <+23219>:	movaps 0xfa0(%r8),%xmm15
   0x0000000000005ad7 <+23255>:	movaps 0x3e0(%rcx),%xmm9
   0x0000000000005adf <+23263>:	mulps  %xmm9,%xmm8
   0x0000000000005ae6 <+23270>:	mulps  %xmm9,%xmm15
   0x0000000000005aea <+23274>:	movaps 0xf60(%r8),%xmm4
   0x0000000000005af2 <+23282>:	addps  %xmm7,%xmm8
   0x0000000000005af6 <+23286>:	mulps  %xmm9,%xmm4
   0x0000000000005b10 <+23312>:	addps  %xmm3,%xmm4
   0x0000000000005b4e <+23374>:	movaps 0xfe0(%r8),%xmm10
   0x0000000000005b5a <+23386>:	mulps  %xmm9,%xmm10
   0x0000000000005b66 <+23398>:	addps  %xmm5,%xmm15
   0x0000000000005b75 <+23413>:	addps  %xmm3,%xmm10

1757
1758	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1759	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
1760	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005a5d <+23133>:	movaps 0xf30(%r8),%xmm12
   0x0000000000005abb <+23227>:	movaps 0xfb0(%r8),%xmm6
   0x0000000000005afa <+23290>:	movaps 0x3f0(%rcx),%xmm7
   0x0000000000005b04 <+23300>:	mulps  %xmm7,%xmm12
   0x0000000000005b1b <+23323>:	mulps  %xmm7,%xmm6
   0x0000000000005b1e <+23326>:	addps  %xmm8,%xmm12
   0x0000000000005b26 <+23334>:	movaps 0xf70(%r8),%xmm8
   0x0000000000005b3a <+23354>:	mulps  %xmm7,%xmm8
   0x0000000000005b46 <+23366>:	addps  %xmm4,%xmm8
   0x0000000000005b5e <+23390>:	movaps 0xff0(%r8),%xmm9
   0x0000000000005b6a <+23402>:	mulps  %xmm7,%xmm9
   0x0000000000005b6e <+23406>:	addps  %xmm15,%xmm6
   0x0000000000005b79 <+23417>:	addps  %xmm10,%xmm9

1761	      }
1762	    }
1763
1764	    /* Store C(4,4) block. */
1765	    for (i = 0; i < 4; i++)
1766	    {
1767	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000005b7d <+23421>:	mulps  %xmm0,%xmm12
   0x0000000000005b99 <+23449>:	mulps  %xmm0,%xmm8
   0x0000000000005bad <+23469>:	mulps  %xmm0,%xmm6
   0x0000000000005bc0 <+23488>:	mulps  %xmm0,%xmm9

1768	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_44]), C_row[i]);
   0x0000000000005b81 <+23425>:	addps  0x3c0(%r9),%xmm12
   0x0000000000005b9d <+23453>:	addps  0x3d0(%r9),%xmm8
   0x0000000000005bb0 <+23472>:	addps  0x3e0(%r9),%xmm6
   0x0000000000005bc4 <+23492>:	addps  0x3f0(%r9),%xmm9

1769	      _mm_store_ps(&C[i*4+C_OFFSET_44], C_row[i]);
   0x0000000000005b8d <+23437>:	movaps %xmm12,0x3c0(%r9)
   0x0000000000005ba5 <+23461>:	movaps %xmm8,0x3d0(%r9)
   0x0000000000005bb8 <+23480>:	movaps %xmm6,0x3e0(%r9)
   0x0000000000005bcc <+23500>:	movaps %xmm9,0x3f0(%r9)

1770	    }
1771	  }
1772	}
   0x0000000000005bda <+23514>:	add    $0x1e8,%rsp
   0x0000000000005be1 <+23521>:	retq
   0x0000000000005be2 <+23522>:	nopl   0x0(%rax)
   0x0000000000005be9 <+23529>:	nopl   0x0(%rax)

End of assembler dump.
