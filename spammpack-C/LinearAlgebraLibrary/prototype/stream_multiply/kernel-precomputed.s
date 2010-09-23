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
   0x0000000000000000 <+0>:	sub    $0x208,%rsp

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
   0x0000000000000012 <+18>:	jbe    0x5ac3 <stream_kernel+23235>
   0x0000000000000018 <+24>:	xor    %eax,%eax
   0x000000000000001a <+26>:	movaps %xmm0,0x20(%rsp)
   0x000000000000001f <+31>:	movss  %xmm1,0x1b8(%rsp)
   0x0000000000005ab5 <+23221>:	inc    %eax
   0x0000000000005ab7 <+23223>:	mov    %eax,%edx
   0x0000000000005ab9 <+23225>:	mov    %edx,%eax
   0x0000000000005abb <+23227>:	cmp    %edi,%eax
   0x0000000000005abd <+23229>:	jb     0x1f <stream_kernel+31>

99	  {
100	    /* Load pointers to matrix data blocks. */
101	    A = multiply_stream[stream_index].A_block;
   0x0000000000000028 <+40>:	imul   $0x98,%rdx,%rdx
   0x000000000000002f <+47>:	mov    (%rsi,%rdx,1),%r8

102	    B = multiply_stream[stream_index].B_block;
   0x0000000000000033 <+51>:	mov    0x8(%rsi,%rdx,1),%rcx

103	    C = multiply_stream[stream_index].C_block;
   0x0000000000000038 <+56>:	mov    0x10(%rsi,%rdx,1),%r9

104	    norm = multiply_stream[stream_index].norm;
105
106	    /* Calculate norms. */
107	    norm_product[0][0] = (norm[0]*norm[16] >= tolerance);
   0x000000000000003d <+61>:	movss  0x18(%rdx,%rsi,1),%xmm3
   0x0000000000000043 <+67>:	movss  0x58(%rdx,%rsi,1),%xmm1
   0x000000000000005c <+92>:	movaps %xmm3,%xmm6
   0x000000000000005f <+95>:	mulss  %xmm1,%xmm6
   0x0000000000000081 <+129>:	movss  %xmm6,0x1c0(%rsp)

108	    norm_product[1][4] = (norm[1]*norm[20] >= tolerance);
   0x0000000000000049 <+73>:	movss  0x1c(%rdx,%rsi,1),%xmm2
   0x0000000000000063 <+99>:	movss  0x68(%rdx,%rsi,1),%xmm13
   0x000000000000008a <+138>:	movaps %xmm2,%xmm5
   0x000000000000008d <+141>:	mulss  %xmm13,%xmm5
   0x000000000000009d <+157>:	movss  %xmm5,0x1b0(%rsp)

109	    norm_product[2][8] = (norm[2]*norm[24] >= tolerance);
   0x000000000000004f <+79>:	movss  0x20(%rdx,%rsi,1),%xmm0
   0x000000000000006a <+106>:	movss  0x78(%rdx,%rsi,1),%xmm7
   0x0000000000000092 <+146>:	movaps %xmm0,%xmm4
   0x0000000000000099 <+153>:	mulss  %xmm7,%xmm4
   0x00000000000000b6 <+182>:	movss  %xmm4,0x1a8(%rsp)

110	    norm_product[3][12] = (norm[3]*norm[28] >= tolerance);
   0x0000000000000055 <+85>:	movss  0x24(%rdx,%rsi,1),%xmm15
   0x0000000000000070 <+112>:	movss  0x88(%rdx,%rsi,1),%xmm8
   0x0000000000000095 <+149>:	movaps %xmm15,%xmm6
   0x00000000000000b1 <+177>:	mulss  %xmm8,%xmm6
   0x00000000000000cf <+207>:	movss  %xmm6,0x1a0(%rsp)

111	    norm_product[0][1] = (norm[0]*norm[17] >= tolerance);
   0x00000000000000a6 <+166>:	movaps %xmm3,%xmm5
   0x00000000000000d8 <+216>:	movss  0x5c(%rdx,%rsi,1),%xmm6
   0x00000000000000de <+222>:	mulss  %xmm6,%xmm5
   0x00000000000000f8 <+248>:	movss  %xmm5,0x198(%rsp)

112	    norm_product[1][5] = (norm[1]*norm[21] >= tolerance);
   0x000000000000007a <+122>:	movss  0x6c(%rdx,%rsi,1),%xmm9
   0x00000000000000a9 <+169>:	movaps %xmm2,%xmm12
   0x00000000000000ca <+202>:	mulss  %xmm9,%xmm12
   0x00000000000000e2 <+226>:	movss  %xmm12,0x190(%rsp)

113	    norm_product[2][9] = (norm[2]*norm[25] >= tolerance);
   0x00000000000000ad <+173>:	movaps %xmm0,%xmm14
   0x00000000000000ec <+236>:	movss  0x7c(%rdx,%rsi,1),%xmm12
   0x00000000000000f3 <+243>:	mulss  %xmm12,%xmm14
   0x000000000000010c <+268>:	movss  %xmm14,0x188(%rsp)

114	    norm_product[3][13] = (norm[3]*norm[29] >= tolerance);
   0x00000000000000bf <+191>:	movaps %xmm15,%xmm11
   0x0000000000000116 <+278>:	movss  0x8c(%rdx,%rsi,1),%xmm14
   0x0000000000000120 <+288>:	mulss  %xmm14,%xmm11
   0x0000000000000139 <+313>:	movss  %xmm11,0x180(%rsp)

115	    norm_product[0][2] = (norm[0]*norm[18] >= tolerance);
   0x00000000000000c3 <+195>:	movaps %xmm3,%xmm10
   0x0000000000000101 <+257>:	movss  0x60(%rdx,%rsi,1),%xmm5
   0x0000000000000107 <+263>:	mulss  %xmm5,%xmm10
   0x0000000000000125 <+293>:	movss  %xmm10,0x1c8(%rsp)

116	    norm_product[1][6] = (norm[1]*norm[22] >= tolerance);
   0x00000000000000c7 <+199>:	movaps %xmm2,%xmm4
   0x0000000000000143 <+323>:	movss  0x70(%rdx,%rsi,1),%xmm11
   0x000000000000014a <+330>:	mulss  %xmm11,%xmm4
   0x000000000000014f <+335>:	movss  %xmm4,0x170(%rsp)

117	    norm_product[2][10] = (norm[2]*norm[26] >= tolerance);
   0x000000000000012f <+303>:	movss  0x80(%rdx,%rsi,1),%xmm10
   0x0000000000000158 <+344>:	movaps %xmm0,%xmm4
   0x000000000000015b <+347>:	mulss  %xmm10,%xmm4
   0x0000000000000160 <+352>:	movss  %xmm4,0x168(%rsp)

118	    norm_product[3][14] = (norm[3]*norm[30] >= tolerance);
   0x0000000000000169 <+361>:	movaps %xmm15,%xmm4
   0x000000000000016d <+365>:	mulss  0x90(%rdx,%rsi,1),%xmm4
   0x0000000000000176 <+374>:	movss  %xmm4,0x158(%rsp)

119	    norm_product[0][3] = (norm[0]*norm[19] >= tolerance);
   0x000000000000017f <+383>:	movss  0x64(%rdx,%rsi,1),%xmm4
   0x0000000000000185 <+389>:	mulss  %xmm4,%xmm3
   0x0000000000000189 <+393>:	movss  %xmm3,0x1d0(%rsp)

120	    norm_product[1][7] = (norm[1]*norm[23] >= tolerance);
   0x0000000000000192 <+402>:	movss  0x74(%rdx,%rsi,1),%xmm3
   0x0000000000000198 <+408>:	mulss  %xmm3,%xmm2
   0x000000000000019c <+412>:	movss  %xmm2,0x150(%rsp)

121	    norm_product[2][11] = (norm[2]*norm[27] >= tolerance);
   0x00000000000001a5 <+421>:	movss  0x84(%rdx,%rsi,1),%xmm2
   0x00000000000001ae <+430>:	mulss  %xmm2,%xmm0
   0x00000000000001b2 <+434>:	movss  %xmm0,0x140(%rsp)

122	    norm_product[3][15] = (norm[3]*norm[31] >= tolerance);
   0x00000000000001bb <+443>:	movss  0x94(%rdx,%rsi,1),%xmm0
   0x00000000000001c4 <+452>:	mulss  %xmm0,%xmm15
   0x00000000000001c9 <+457>:	movss  %xmm15,0x138(%rsp)

123	    norm_product[4][0] = (norm[4]*norm[16] >= tolerance);
   0x00000000000001d3 <+467>:	movss  0x28(%rdx,%rsi,1),%xmm15
   0x00000000000001da <+474>:	mulss  %xmm15,%xmm1
   0x00000000000001df <+479>:	movss  %xmm1,0x148(%rsp)

124	    norm_product[5][4] = (norm[5]*norm[20] >= tolerance);
   0x00000000000001e8 <+488>:	movss  0x2c(%rdx,%rsi,1),%xmm1
   0x00000000000001ee <+494>:	mulss  %xmm1,%xmm13
   0x00000000000001f3 <+499>:	movss  %xmm13,0x130(%rsp)

125	    norm_product[6][8] = (norm[6]*norm[24] >= tolerance);
   0x00000000000001fd <+509>:	movaps %xmm7,%xmm13
   0x0000000000000201 <+513>:	movss  0x30(%rdx,%rsi,1),%xmm7
   0x0000000000000207 <+519>:	mulss  %xmm7,%xmm13
   0x000000000000020c <+524>:	movss  %xmm13,0x120(%rsp)

126	    norm_product[7][12] = (norm[7]*norm[28] >= tolerance);
   0x0000000000000216 <+534>:	movss  0x34(%rdx,%rsi,1),%xmm13
   0x000000000000021d <+541>:	mulss  %xmm13,%xmm8
   0x0000000000000222 <+546>:	movss  %xmm8,0x118(%rsp)

127	    norm_product[4][1] = (norm[4]*norm[17] >= tolerance);
   0x000000000000022c <+556>:	movaps %xmm6,%xmm8
   0x0000000000000230 <+560>:	mulss  %xmm15,%xmm8
   0x0000000000000235 <+565>:	movss  %xmm8,0x1d8(%rsp)

128	    norm_product[5][5] = (norm[5]*norm[21] >= tolerance);
   0x000000000000023f <+575>:	movaps %xmm9,%xmm8
   0x0000000000000243 <+579>:	mulss  %xmm1,%xmm8
   0x0000000000000248 <+584>:	movss  %xmm8,0x110(%rsp)

129	    norm_product[6][9] = (norm[6]*norm[25] >= tolerance);
   0x0000000000000252 <+594>:	movaps %xmm12,%xmm8
   0x0000000000000256 <+598>:	mulss  %xmm7,%xmm8
   0x000000000000025b <+603>:	movss  %xmm8,0x108(%rsp)

130	    norm_product[7][13] = (norm[7]*norm[29] >= tolerance);
   0x0000000000000265 <+613>:	movaps %xmm14,%xmm8
   0x0000000000000269 <+617>:	mulss  %xmm13,%xmm8
   0x000000000000026e <+622>:	movss  %xmm8,0x100(%rsp)

131	    norm_product[4][2] = (norm[4]*norm[18] >= tolerance);
   0x0000000000000278 <+632>:	movaps %xmm5,%xmm8
   0x000000000000027c <+636>:	mulss  %xmm15,%xmm8
   0x0000000000000286 <+646>:	movss  %xmm8,0x178(%rsp)

132	    norm_product[5][6] = (norm[5]*norm[22] >= tolerance);
   0x0000000000000290 <+656>:	movaps %xmm11,%xmm8
   0x0000000000000294 <+660>:	mulss  %xmm1,%xmm8
   0x00000000000002ae <+686>:	movss  %xmm8,0xf8(%rsp)

133	    norm_product[6][10] = (norm[6]*norm[26] >= tolerance);
   0x00000000000002b8 <+696>:	movaps %xmm10,%xmm8
   0x00000000000002bc <+700>:	mulss  %xmm7,%xmm8
   0x00000000000002d4 <+724>:	movss  %xmm8,0xe8(%rsp)

134	    norm_product[7][14] = (norm[7]*norm[30] >= tolerance);
   0x00000000000002de <+734>:	movss  0x90(%rdx,%rsi,1),%xmm8
   0x00000000000002e8 <+744>:	mulss  %xmm13,%xmm8
   0x0000000000000301 <+769>:	movss  %xmm8,0xe0(%rsp)

135	    norm_product[4][3] = (norm[4]*norm[19] >= tolerance);
   0x0000000000000281 <+641>:	mulss  %xmm4,%xmm15
   0x0000000000000299 <+665>:	movss  %xmm15,0x1e0(%rsp)

136	    norm_product[5][7] = (norm[5]*norm[23] >= tolerance);
   0x00000000000002aa <+682>:	mulss  %xmm3,%xmm1
   0x00000000000002c1 <+705>:	movss  %xmm1,0xf0(%rsp)

137	    norm_product[6][11] = (norm[6]*norm[27] >= tolerance);
   0x00000000000002d0 <+720>:	mulss  %xmm2,%xmm7
   0x00000000000002ed <+749>:	movss  %xmm7,0xd8(%rsp)

138	    norm_product[7][15] = (norm[7]*norm[31] >= tolerance);
   0x00000000000002fc <+764>:	mulss  %xmm0,%xmm13
   0x0000000000000317 <+791>:	movss  %xmm13,0xd0(%rsp)

139	    norm_product[8][0] = (norm[8]*norm[16] >= tolerance);
   0x00000000000002a3 <+675>:	movss  0x58(%rdx,%rsi,1),%xmm15
   0x00000000000002ca <+714>:	movss  0x38(%rdx,%rsi,1),%xmm1
   0x0000000000000312 <+786>:	mulss  %xmm1,%xmm15
   0x000000000000032d <+813>:	movss  %xmm15,0x160(%rsp)

140	    norm_product[9][4] = (norm[9]*norm[20] >= tolerance);
   0x0000000000000321 <+801>:	movss  0x68(%rdx,%rsi,1),%xmm13
   0x0000000000000337 <+823>:	movss  0x3c(%rdx,%rsi,1),%xmm15
   0x000000000000033e <+830>:	mulss  %xmm15,%xmm13
   0x0000000000000352 <+850>:	movss  %xmm13,0xc8(%rsp)

141	    norm_product[10][8] = (norm[10]*norm[24] >= tolerance);
   0x00000000000002f6 <+758>:	movss  0x78(%rdx,%rsi,1),%xmm7
   0x000000000000030b <+779>:	movss  0x40(%rdx,%rsi,1),%xmm8
   0x0000000000000328 <+808>:	mulss  %xmm8,%xmm7
   0x0000000000000343 <+835>:	movss  %xmm7,0xc0(%rsp)

142	    norm_product[11][12] = (norm[11]*norm[28] >= tolerance);
   0x000000000000034c <+844>:	movss  0x44(%rdx,%rsi,1),%xmm7
   0x000000000000035c <+860>:	movss  0x88(%rdx,%rsi,1),%xmm13
   0x0000000000000366 <+870>:	mulss  %xmm7,%xmm13
   0x000000000000036b <+875>:	movss  %xmm13,0xa8(%rsp)

143	    norm_product[8][1] = (norm[8]*norm[17] >= tolerance);
   0x0000000000000375 <+885>:	movaps %xmm6,%xmm13
   0x0000000000000379 <+889>:	mulss  %xmm1,%xmm13
   0x000000000000037e <+894>:	movss  %xmm13,0x128(%rsp)

144	    norm_product[9][5] = (norm[9]*norm[21] >= tolerance);
   0x0000000000000388 <+904>:	movaps %xmm9,%xmm13
   0x000000000000038c <+908>:	mulss  %xmm15,%xmm13
   0x0000000000000391 <+913>:	movss  %xmm13,0xb0(%rsp)

145	    norm_product[10][9] = (norm[10]*norm[25] >= tolerance);
   0x000000000000039b <+923>:	movaps %xmm12,%xmm13
   0x000000000000039f <+927>:	mulss  %xmm8,%xmm13
   0x00000000000003a4 <+932>:	movss  %xmm13,0x98(%rsp)

146	    norm_product[11][13] = (norm[11]*norm[29] >= tolerance);
   0x00000000000003ae <+942>:	movaps %xmm14,%xmm13
   0x00000000000003b2 <+946>:	mulss  %xmm7,%xmm13
   0x00000000000003b7 <+951>:	movss  %xmm13,0x88(%rsp)

147	    norm_product[8][2] = (norm[8]*norm[18] >= tolerance);
   0x00000000000003c1 <+961>:	movaps %xmm5,%xmm13
   0x00000000000003c5 <+965>:	mulss  %xmm1,%xmm13
   0x00000000000003ce <+974>:	movss  %xmm13,0x1e8(%rsp)

148	    norm_product[9][6] = (norm[9]*norm[22] >= tolerance);
   0x00000000000003d8 <+984>:	movaps %xmm11,%xmm13
   0x00000000000003dc <+988>:	mulss  %xmm15,%xmm13
   0x00000000000003f5 <+1013>:	movss  %xmm13,0x90(%rsp)

149	    norm_product[10][10] = (norm[10]*norm[26] >= tolerance);
   0x00000000000003ff <+1023>:	movaps %xmm10,%xmm13
   0x0000000000000403 <+1027>:	mulss  %xmm8,%xmm13
   0x000000000000041b <+1051>:	movss  %xmm13,0x70(%rsp)

150	    norm_product[11][14] = (norm[11]*norm[30] >= tolerance);
   0x0000000000000422 <+1058>:	movss  0x90(%rdx,%rsi,1),%xmm13
   0x000000000000042c <+1068>:	mulss  %xmm7,%xmm13
   0x0000000000000443 <+1091>:	movss  %xmm13,0x60(%rsp)

151	    norm_product[8][3] = (norm[8]*norm[19] >= tolerance);
   0x00000000000003ca <+970>:	mulss  %xmm4,%xmm1
   0x00000000000003e1 <+993>:	movss  %xmm1,0x1f0(%rsp)

152	    norm_product[9][7] = (norm[9]*norm[23] >= tolerance);
   0x00000000000003f0 <+1008>:	mulss  %xmm3,%xmm15
   0x0000000000000408 <+1032>:	movss  %xmm15,0x78(%rsp)

153	    norm_product[10][11] = (norm[10]*norm[27] >= tolerance);
   0x0000000000000416 <+1046>:	mulss  %xmm2,%xmm8
   0x0000000000000431 <+1073>:	movss  %xmm8,0x68(%rsp)

154	    norm_product[11][15] = (norm[11]*norm[31] >= tolerance);
   0x000000000000043f <+1087>:	mulss  %xmm0,%xmm7
   0x0000000000000456 <+1110>:	movss  %xmm7,0x40(%rsp)

155	    norm_product[12][0] = (norm[12]*norm[16] >= tolerance);
   0x00000000000003ea <+1002>:	movss  0x48(%rdx,%rsi,1),%xmm1
   0x0000000000000438 <+1080>:	movss  0x58(%rdx,%rsi,1),%xmm8
   0x0000000000000451 <+1105>:	mulss  %xmm1,%xmm8
   0x0000000000000467 <+1127>:	movss  %xmm8,0x1f8(%rsp)

156	    norm_product[13][4] = (norm[13]*norm[20] >= tolerance);
   0x000000000000040f <+1039>:	movss  0x68(%rdx,%rsi,1),%xmm15
   0x000000000000044a <+1098>:	movss  0x4c(%rdx,%rsi,1),%xmm13
   0x0000000000000462 <+1122>:	mulss  %xmm13,%xmm15
   0x000000000000047d <+1149>:	movss  %xmm15,0xb8(%rsp)

157	    norm_product[14][8] = (norm[14]*norm[24] >= tolerance);
   0x000000000000045c <+1116>:	movss  0x50(%rdx,%rsi,1),%xmm7
   0x0000000000000471 <+1137>:	movss  0x78(%rdx,%rsi,1),%xmm8
   0x0000000000000478 <+1144>:	mulss  %xmm7,%xmm8
   0x0000000000000495 <+1173>:	movss  %xmm8,0x58(%rsp)

158	    norm_product[15][12] = (norm[15]*norm[28] >= tolerance);
   0x0000000000000487 <+1159>:	movss  0x88(%rdx,%rsi,1),%xmm15
   0x000000000000049c <+1180>:	movss  0x54(%rdx,%rsi,1),%xmm8
   0x00000000000004a3 <+1187>:	mulss  %xmm8,%xmm15
   0x00000000000004a8 <+1192>:	movss  %xmm15,0x38(%rsp)

159	    norm_product[12][1] = (norm[12]*norm[17] >= tolerance);
   0x0000000000000491 <+1169>:	mulss  %xmm1,%xmm6

160	    norm_product[13][5] = (norm[13]*norm[21] >= tolerance);
   0x00000000000004af <+1199>:	mulss  %xmm13,%xmm9
   0x00000000000004b4 <+1204>:	movss  %xmm9,0xa0(%rsp)

161	    norm_product[14][9] = (norm[14]*norm[25] >= tolerance);
   0x00000000000004c8 <+1224>:	mulss  %xmm7,%xmm12
   0x00000000000004cd <+1229>:	movss  %xmm12,0x50(%rsp)

162	    norm_product[15][13] = (norm[15]*norm[29] >= tolerance);
   0x00000000000004de <+1246>:	mulss  %xmm8,%xmm14
   0x00000000000004e3 <+1251>:	movss  %xmm14,0x18(%rsp)

163	    norm_product[12][2] = (norm[12]*norm[18] >= tolerance);
   0x00000000000004ea <+1258>:	mulss  %xmm1,%xmm5

164	    norm_product[13][6] = (norm[13]*norm[22] >= tolerance);
   0x00000000000004ee <+1262>:	mulss  %xmm13,%xmm11
   0x00000000000004f3 <+1267>:	movss  %xmm11,0x80(%rsp)

165	    norm_product[14][10] = (norm[14]*norm[26] >= tolerance);
   0x0000000000000507 <+1287>:	mulss  %xmm7,%xmm10
   0x000000000000050c <+1292>:	movss  %xmm10,0x48(%rsp)

166	    norm_product[15][14] = (norm[15]*norm[30] >= tolerance);
   0x00000000000004be <+1214>:	movss  0x90(%rdx,%rsi,1),%xmm9
   0x000000000000051d <+1309>:	mulss  %xmm8,%xmm9
   0x0000000000000522 <+1314>:	movss  %xmm9,0x30(%rsp)

167	    norm_product[12][3] = (norm[12]*norm[19] >= tolerance);
   0x0000000000000533 <+1331>:	mulss  %xmm1,%xmm4

168	    norm_product[13][7] = (norm[13]*norm[23] >= tolerance);
   0x0000000000000540 <+1344>:	mulss  %xmm13,%xmm3

169	    norm_product[14][11] = (norm[14]*norm[27] >= tolerance);
   0x000000000000054f <+1359>:	mulss  %xmm7,%xmm2

170	    norm_product[15][15] = (norm[15]*norm[31] >= tolerance);
   0x0000000000000568 <+1384>:	mulss  %xmm8,%xmm0

171
172	    /* Reset C(1,1) matrix accumulators */
173	    C_row[0] = _mm_setzero_ps();
174	    C_row[1] = _mm_setzero_ps();
175	    C_row[2] = _mm_setzero_ps();
176	    C_row[3] = _mm_setzero_ps();
177
178	    if (norm_product[0][0] &&
   0x00000000000004d4 <+1236>:	movss  0x1d0(%rsp),%xmm12
   0x00000000000004fd <+1277>:	movss  0x1d8(%rsp),%xmm11
   0x0000000000000513 <+1299>:	movss  0x1e0(%rsp),%xmm10
   0x0000000000000529 <+1321>:	movss  0x1e8(%rsp),%xmm9
   0x0000000000000537 <+1335>:	movss  0x1b8(%rsp),%xmm1
   0x0000000000000545 <+1349>:	movss  0x1c8(%rsp),%xmm13
   0x0000000000000553 <+1363>:	movss  0x1c0(%rsp),%xmm7
   0x000000000000055c <+1372>:	comiss %xmm1,%xmm7
   0x000000000000055f <+1375>:	movss  0x1f8(%rsp),%xmm7
   0x000000000000056d <+1389>:	movss  0x1f0(%rsp),%xmm8
   0x0000000000000577 <+1399>:	jb     0xa9e <stream_kernel+2718>

179	        norm_product[1][4] &&
   0x000000000000057d <+1405>:	movss  0x1b0(%rsp),%xmm14
   0x0000000000000587 <+1415>:	comiss %xmm1,%xmm14
   0x000000000000058b <+1419>:	jb     0xa9e <stream_kernel+2718>

180	        norm_product[2][8] &&
   0x0000000000000591 <+1425>:	movss  0x1a8(%rsp),%xmm14
   0x000000000000059b <+1435>:	comiss %xmm1,%xmm14
   0x000000000000059f <+1439>:	jb     0xa9e <stream_kernel+2718>

181	        norm_product[3][12])
   0x00000000000005a5 <+1445>:	movss  0x1a0(%rsp),%xmm14
   0x00000000000005af <+1455>:	comiss %xmm1,%xmm14
   0x00000000000005b3 <+1459>:	jb     0xa9e <stream_kernel+2718>

182	    {
183	      /* A(1,1)*B(1,1) = C(1,1). */
184	      for (i = 0; i < 4; i++)
185	      {
186	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
187	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
188	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000005b9 <+1465>:	movaps (%r8),%xmm15
   0x00000000000005bd <+1469>:	mulps  (%rcx),%xmm15
   0x00000000000005df <+1503>:	movaps 0x40(%r8),%xmm2
   0x00000000000005e4 <+1508>:	mulps  (%rcx),%xmm2
   0x00000000000005e7 <+1511>:	movss  %xmm0,0x10(%rsp)
   0x0000000000000608 <+1544>:	movaps 0x80(%r8),%xmm2
   0x0000000000000610 <+1552>:	mulps  (%rcx),%xmm2
   0x0000000000000672 <+1650>:	movaps 0xc0(%r8),%xmm14
   0x000000000000067a <+1658>:	mulps  (%rcx),%xmm14

189
190	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
191	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
192	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000005c1 <+1473>:	movaps 0x50(%r8),%xmm14
   0x00000000000005c6 <+1478>:	mulps  0x10(%rcx),%xmm14
   0x00000000000005cb <+1483>:	movss  %xmm3,(%rsp)
   0x00000000000005d0 <+1488>:	movaps 0x10(%r8),%xmm3
   0x00000000000005d5 <+1493>:	mulps  0x10(%rcx),%xmm3
   0x00000000000005d9 <+1497>:	movss  %xmm2,0x8(%rsp)
   0x00000000000005ed <+1517>:	addps  %xmm15,%xmm3
   0x0000000000000604 <+1540>:	addps  %xmm2,%xmm14
   0x0000000000000643 <+1603>:	movaps 0x90(%r8),%xmm3
   0x000000000000064b <+1611>:	mulps  0x10(%rcx),%xmm3
   0x000000000000064f <+1615>:	addps  %xmm2,%xmm3
   0x0000000000000652 <+1618>:	movaps 0xd0(%r8),%xmm2
   0x000000000000065a <+1626>:	mulps  0x10(%rcx),%xmm2
   0x000000000000067e <+1662>:	addps  %xmm14,%xmm2

193
194	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
195	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
196	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000005f1 <+1521>:	movaps 0x20(%r8),%xmm15
   0x00000000000005fb <+1531>:	mulps  0x20(%rcx),%xmm15
   0x0000000000000613 <+1555>:	addps  %xmm3,%xmm15
   0x0000000000000617 <+1559>:	movaps 0x60(%r8),%xmm3
   0x000000000000061c <+1564>:	mulps  0x20(%rcx),%xmm3
   0x000000000000062e <+1582>:	addps  %xmm14,%xmm3
   0x0000000000000632 <+1586>:	movaps 0xa0(%r8),%xmm14
   0x000000000000063a <+1594>:	mulps  0x20(%rcx),%xmm14
   0x000000000000065e <+1630>:	addps  %xmm3,%xmm14
   0x0000000000000682 <+1666>:	movaps 0xe0(%r8),%xmm14
   0x000000000000068a <+1674>:	mulps  0x20(%rcx),%xmm14
   0x000000000000068f <+1679>:	addps  %xmm2,%xmm14

197
198	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
199	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
200	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000005f6 <+1526>:	movaps 0x30(%r8),%xmm0
   0x0000000000000600 <+1536>:	mulps  0x30(%rcx),%xmm0
   0x0000000000000620 <+1568>:	addps  %xmm15,%xmm0
   0x0000000000000624 <+1572>:	movaps 0x70(%r8),%xmm15
   0x0000000000000629 <+1577>:	mulps  0x30(%rcx),%xmm15
   0x000000000000063f <+1599>:	addps  %xmm3,%xmm15
   0x0000000000000662 <+1634>:	movaps 0xb0(%r8),%xmm3
   0x000000000000066a <+1642>:	mulps  0x30(%rcx),%xmm3
   0x000000000000066e <+1646>:	addps  %xmm14,%xmm3
   0x0000000000000693 <+1683>:	movaps 0xf0(%r8),%xmm2
   0x000000000000069b <+1691>:	mulps  0x30(%rcx),%xmm2
   0x000000000000069f <+1695>:	addps  %xmm14,%xmm2

201	      }
202
203	      /* A(1,2)*B(2,1) = C(1,1). */
204	      for (i = 0; i < 4; i++)
205	      {
206	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
207	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
208	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000006a3 <+1699>:	movaps 0x100(%r8),%xmm14
   0x00000000000006ab <+1707>:	mulps  0x100(%rcx),%xmm14
   0x00000000000006b3 <+1715>:	addps  %xmm0,%xmm14
   0x00000000000006f1 <+1777>:	movaps 0x140(%r8),%xmm14
   0x00000000000006f9 <+1785>:	mulps  0x100(%rcx),%xmm14
   0x0000000000000701 <+1793>:	addps  %xmm15,%xmm14
   0x0000000000000741 <+1857>:	movaps 0x180(%r8),%xmm14
   0x0000000000000749 <+1865>:	mulps  0x100(%rcx),%xmm14
   0x0000000000000751 <+1873>:	addps  %xmm3,%xmm14
   0x000000000000078f <+1935>:	movaps 0x1c0(%r8),%xmm14
   0x0000000000000797 <+1943>:	mulps  0x100(%rcx),%xmm14
   0x000000000000079f <+1951>:	addps  %xmm2,%xmm14

209
210	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
211	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
212	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000006b7 <+1719>:	movaps 0x110(%r8),%xmm0
   0x00000000000006bf <+1727>:	mulps  0x110(%rcx),%xmm0
   0x00000000000006c6 <+1734>:	addps  %xmm14,%xmm0
   0x0000000000000705 <+1797>:	movaps 0x150(%r8),%xmm15
   0x000000000000070d <+1805>:	mulps  0x110(%rcx),%xmm15
   0x0000000000000715 <+1813>:	addps  %xmm14,%xmm15
   0x0000000000000755 <+1877>:	movaps 0x190(%r8),%xmm3
   0x000000000000075d <+1885>:	mulps  0x110(%rcx),%xmm3
   0x0000000000000764 <+1892>:	addps  %xmm14,%xmm3
   0x00000000000007a3 <+1955>:	movaps 0x1d0(%r8),%xmm2
   0x00000000000007ab <+1963>:	mulps  0x110(%rcx),%xmm2
   0x00000000000007b2 <+1970>:	addps  %xmm14,%xmm2

213
214	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
215	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
216	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000006ca <+1738>:	movaps 0x120(%r8),%xmm14
   0x00000000000006d2 <+1746>:	mulps  0x120(%rcx),%xmm14
   0x00000000000006da <+1754>:	addps  %xmm0,%xmm14
   0x0000000000000719 <+1817>:	movaps 0x160(%r8),%xmm14
   0x0000000000000721 <+1825>:	mulps  0x120(%rcx),%xmm14
   0x0000000000000729 <+1833>:	addps  %xmm15,%xmm14
   0x0000000000000768 <+1896>:	movaps 0x1a0(%r8),%xmm14
   0x0000000000000770 <+1904>:	mulps  0x120(%rcx),%xmm14
   0x0000000000000778 <+1912>:	addps  %xmm3,%xmm14
   0x00000000000007b6 <+1974>:	movaps 0x1e0(%r8),%xmm14
   0x00000000000007be <+1982>:	mulps  0x120(%rcx),%xmm14
   0x00000000000007c6 <+1990>:	addps  %xmm2,%xmm14

217
218	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
219	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
220	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000006de <+1758>:	movaps 0x130(%r8),%xmm0
   0x00000000000006e6 <+1766>:	mulps  0x130(%rcx),%xmm0
   0x00000000000006ed <+1773>:	addps  %xmm14,%xmm0
   0x000000000000072d <+1837>:	movaps 0x170(%r8),%xmm15
   0x0000000000000735 <+1845>:	mulps  0x130(%rcx),%xmm15
   0x000000000000073d <+1853>:	addps  %xmm14,%xmm15
   0x000000000000077c <+1916>:	movaps 0x1b0(%r8),%xmm3
   0x0000000000000784 <+1924>:	mulps  0x130(%rcx),%xmm3
   0x000000000000078b <+1931>:	addps  %xmm14,%xmm3
   0x00000000000007ca <+1994>:	movaps 0x1f0(%r8),%xmm2
   0x00000000000007d2 <+2002>:	mulps  0x130(%rcx),%xmm2
   0x00000000000007d9 <+2009>:	addps  %xmm14,%xmm2

221	      }
222
223	      /* A(1,3)*B(3,1) = C(1,1). */
224	      for (i = 0; i < 4; i++)
225	      {
226	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
227	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
228	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000007dd <+2013>:	movaps 0x200(%r8),%xmm14
   0x00000000000007e5 <+2021>:	mulps  0x200(%rcx),%xmm14
   0x00000000000007ed <+2029>:	addps  %xmm0,%xmm14
   0x000000000000082b <+2091>:	movaps 0x240(%r8),%xmm14
   0x0000000000000833 <+2099>:	mulps  0x200(%rcx),%xmm14
   0x000000000000083b <+2107>:	addps  %xmm15,%xmm14
   0x000000000000087b <+2171>:	movaps 0x280(%r8),%xmm14
   0x0000000000000883 <+2179>:	mulps  0x200(%rcx),%xmm14
   0x000000000000088b <+2187>:	addps  %xmm3,%xmm14
   0x00000000000008c9 <+2249>:	movaps 0x2c0(%r8),%xmm14
   0x00000000000008d1 <+2257>:	mulps  0x200(%rcx),%xmm14
   0x00000000000008d9 <+2265>:	addps  %xmm2,%xmm14

229
230	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
231	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
232	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000007f1 <+2033>:	movaps 0x210(%r8),%xmm0
   0x00000000000007f9 <+2041>:	mulps  0x210(%rcx),%xmm0
   0x0000000000000800 <+2048>:	addps  %xmm14,%xmm0
   0x000000000000083f <+2111>:	movaps 0x250(%r8),%xmm15
   0x0000000000000847 <+2119>:	mulps  0x210(%rcx),%xmm15
   0x000000000000084f <+2127>:	addps  %xmm14,%xmm15
   0x000000000000088f <+2191>:	movaps 0x290(%r8),%xmm3
   0x0000000000000897 <+2199>:	mulps  0x210(%rcx),%xmm3
   0x000000000000089e <+2206>:	addps  %xmm14,%xmm3
   0x00000000000008dd <+2269>:	movaps 0x2d0(%r8),%xmm2
   0x00000000000008e5 <+2277>:	mulps  0x210(%rcx),%xmm2
   0x00000000000008ec <+2284>:	addps  %xmm14,%xmm2

233
234	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
235	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
236	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000804 <+2052>:	movaps 0x220(%r8),%xmm14
   0x000000000000080c <+2060>:	mulps  0x220(%rcx),%xmm14
   0x0000000000000814 <+2068>:	addps  %xmm0,%xmm14
   0x0000000000000853 <+2131>:	movaps 0x260(%r8),%xmm14
   0x000000000000085b <+2139>:	mulps  0x220(%rcx),%xmm14
   0x0000000000000863 <+2147>:	addps  %xmm15,%xmm14
   0x00000000000008a2 <+2210>:	movaps 0x2a0(%r8),%xmm14
   0x00000000000008aa <+2218>:	mulps  0x220(%rcx),%xmm14
   0x00000000000008b2 <+2226>:	addps  %xmm3,%xmm14
   0x00000000000008f0 <+2288>:	movaps 0x2e0(%r8),%xmm14
   0x00000000000008f8 <+2296>:	mulps  0x220(%rcx),%xmm14
   0x0000000000000900 <+2304>:	addps  %xmm2,%xmm14

237
238	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
239	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
240	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000818 <+2072>:	movaps 0x230(%r8),%xmm0
   0x0000000000000820 <+2080>:	mulps  0x230(%rcx),%xmm0
   0x0000000000000827 <+2087>:	addps  %xmm14,%xmm0
   0x0000000000000867 <+2151>:	movaps 0x270(%r8),%xmm15
   0x000000000000086f <+2159>:	mulps  0x230(%rcx),%xmm15
   0x0000000000000877 <+2167>:	addps  %xmm14,%xmm15
   0x00000000000008b6 <+2230>:	movaps 0x2b0(%r8),%xmm3
   0x00000000000008be <+2238>:	mulps  0x230(%rcx),%xmm3
   0x00000000000008c5 <+2245>:	addps  %xmm14,%xmm3
   0x0000000000000904 <+2308>:	movaps 0x2f0(%r8),%xmm2
   0x000000000000090c <+2316>:	mulps  0x230(%rcx),%xmm2
   0x0000000000000913 <+2323>:	addps  %xmm14,%xmm2

241	      }
242
243	      /* A(1,4)*B(4,1) = C(1,1). */
244	      for (i = 0; i < 4; i++)
245	      {
246	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
247	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
248	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000917 <+2327>:	movaps 0x300(%r8),%xmm14
   0x000000000000091f <+2335>:	mulps  0x300(%rcx),%xmm14
   0x0000000000000927 <+2343>:	addps  %xmm0,%xmm14
   0x0000000000000965 <+2405>:	movaps 0x340(%r8),%xmm14
   0x000000000000096d <+2413>:	mulps  0x300(%rcx),%xmm14
   0x0000000000000975 <+2421>:	addps  %xmm15,%xmm14
   0x00000000000009b5 <+2485>:	movaps 0x380(%r8),%xmm14
   0x00000000000009bd <+2493>:	mulps  0x300(%rcx),%xmm14
   0x00000000000009c5 <+2501>:	addps  %xmm3,%xmm14
   0x0000000000000a03 <+2563>:	movaps 0x3c0(%r8),%xmm14
   0x0000000000000a0b <+2571>:	mulps  0x300(%rcx),%xmm14
   0x0000000000000a13 <+2579>:	addps  %xmm2,%xmm14

249
250	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
251	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
252	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000092b <+2347>:	movaps 0x310(%r8),%xmm0
   0x0000000000000933 <+2355>:	mulps  0x310(%rcx),%xmm0
   0x000000000000093a <+2362>:	addps  %xmm14,%xmm0
   0x0000000000000979 <+2425>:	movaps 0x350(%r8),%xmm15
   0x0000000000000981 <+2433>:	mulps  0x310(%rcx),%xmm15
   0x0000000000000989 <+2441>:	addps  %xmm14,%xmm15
   0x00000000000009c9 <+2505>:	movaps 0x390(%r8),%xmm3
   0x00000000000009d1 <+2513>:	mulps  0x310(%rcx),%xmm3
   0x00000000000009d8 <+2520>:	addps  %xmm14,%xmm3
   0x0000000000000a17 <+2583>:	movaps 0x3d0(%r8),%xmm2
   0x0000000000000a1f <+2591>:	mulps  0x310(%rcx),%xmm2
   0x0000000000000a26 <+2598>:	addps  %xmm14,%xmm2

253
254	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
255	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
256	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000093e <+2366>:	movaps 0x320(%r8),%xmm14
   0x0000000000000946 <+2374>:	mulps  0x320(%rcx),%xmm14
   0x000000000000094e <+2382>:	addps  %xmm0,%xmm14
   0x000000000000098d <+2445>:	movaps 0x360(%r8),%xmm14
   0x0000000000000995 <+2453>:	mulps  0x320(%rcx),%xmm14
   0x000000000000099d <+2461>:	addps  %xmm15,%xmm14
   0x00000000000009dc <+2524>:	movaps 0x3a0(%r8),%xmm14
   0x00000000000009e4 <+2532>:	mulps  0x320(%rcx),%xmm14
   0x00000000000009ec <+2540>:	addps  %xmm3,%xmm14
   0x0000000000000a2a <+2602>:	movaps 0x3e0(%r8),%xmm14
   0x0000000000000a32 <+2610>:	mulps  0x320(%rcx),%xmm14
   0x0000000000000a3a <+2618>:	addps  %xmm2,%xmm14

257
258	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
259	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
260	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000952 <+2386>:	movaps 0x330(%r8),%xmm0
   0x000000000000095a <+2394>:	mulps  0x330(%rcx),%xmm0
   0x0000000000000961 <+2401>:	addps  %xmm14,%xmm0
   0x00000000000009a1 <+2465>:	movaps 0x370(%r8),%xmm15
   0x00000000000009a9 <+2473>:	mulps  0x330(%rcx),%xmm15
   0x00000000000009b1 <+2481>:	addps  %xmm14,%xmm15
   0x00000000000009f0 <+2544>:	movaps 0x3b0(%r8),%xmm3
   0x00000000000009f8 <+2552>:	mulps  0x330(%rcx),%xmm3
   0x00000000000009ff <+2559>:	addps  %xmm14,%xmm3
   0x0000000000000a3e <+2622>:	movaps 0x3f0(%r8),%xmm2
   0x0000000000000a46 <+2630>:	mulps  0x330(%rcx),%xmm2
   0x0000000000000a4d <+2637>:	addps  %xmm14,%xmm2

261	      }
262
263	      /* Store C(1,1) block. */
264	      for (i = 0; i < 4; i++)
265	      {
266	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000000a51 <+2641>:	movaps 0x20(%rsp),%xmm14
   0x0000000000000a57 <+2647>:	mulps  %xmm14,%xmm0
   0x0000000000000a69 <+2665>:	mulps  %xmm14,%xmm15
   0x0000000000000a77 <+2679>:	mulps  %xmm14,%xmm3
   0x0000000000000a8a <+2698>:	mulps  %xmm14,%xmm2

267	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_11]), C_row[i]);
   0x0000000000000a5b <+2651>:	addps  (%r9),%xmm0
   0x0000000000000a6d <+2669>:	addps  0x10(%r9),%xmm15
   0x0000000000000a7b <+2683>:	addps  0x20(%r9),%xmm3
   0x0000000000000a8e <+2702>:	addps  0x30(%r9),%xmm2

268	        _mm_store_ps(&C[i*4+C_OFFSET_11], C_row[i]);
   0x0000000000000a5f <+2655>:	movaps %xmm0,(%r9)
   0x0000000000000a63 <+2659>:	movss  0x10(%rsp),%xmm0
   0x0000000000000a72 <+2674>:	movaps %xmm15,0x10(%r9)
   0x0000000000000a80 <+2688>:	movaps %xmm3,0x20(%r9)
   0x0000000000000a85 <+2693>:	movss  (%rsp),%xmm3
   0x0000000000000a93 <+2707>:	movaps %xmm2,0x30(%r9)
   0x0000000000000a98 <+2712>:	movss  0x8(%rsp),%xmm2

269	      }
270	    }
271
272	    /* Reset C(1,2) matrix accumulators */
273	    C_row[0] = _mm_setzero_ps();
274	    C_row[1] = _mm_setzero_ps();
275	    C_row[2] = _mm_setzero_ps();
276	    C_row[3] = _mm_setzero_ps();
277
278	    if (norm_product[0][1] &&
   0x0000000000000a9e <+2718>:	movss  0x198(%rsp),%xmm14
   0x0000000000000aa8 <+2728>:	comiss %xmm1,%xmm14
   0x0000000000000aac <+2732>:	jb     0xfd9 <stream_kernel+4057>

279	        norm_product[1][5] &&
   0x0000000000000ab2 <+2738>:	movss  0x190(%rsp),%xmm14
   0x0000000000000abc <+2748>:	comiss %xmm1,%xmm14
   0x0000000000000ac0 <+2752>:	jb     0xfd9 <stream_kernel+4057>

280	        norm_product[2][9] &&
   0x0000000000000ac6 <+2758>:	movss  0x188(%rsp),%xmm14
   0x0000000000000ad0 <+2768>:	comiss %xmm1,%xmm14
   0x0000000000000ad4 <+2772>:	jb     0xfd9 <stream_kernel+4057>

281	        norm_product[3][13])
   0x0000000000000ada <+2778>:	movss  0x180(%rsp),%xmm14
   0x0000000000000ae4 <+2788>:	comiss %xmm1,%xmm14
   0x0000000000000ae8 <+2792>:	jb     0xfd9 <stream_kernel+4057>

282	    {
283	      /* A(1,1)*B(1,2) = C(1,2). */
284	      for (i = 0; i < 4; i++)
285	      {
286	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
287	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
288	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000aee <+2798>:	movaps (%r8),%xmm15
   0x0000000000000af2 <+2802>:	mulps  0x40(%rcx),%xmm15
   0x0000000000000b15 <+2837>:	movaps 0x40(%r8),%xmm2
   0x0000000000000b1a <+2842>:	mulps  0x40(%rcx),%xmm2
   0x0000000000000b1e <+2846>:	movss  %xmm0,0x10(%rsp)
   0x0000000000000b3f <+2879>:	movaps 0x80(%r8),%xmm2
   0x0000000000000b47 <+2887>:	mulps  0x40(%rcx),%xmm2
   0x0000000000000baa <+2986>:	movaps 0xc0(%r8),%xmm14
   0x0000000000000bb2 <+2994>:	mulps  0x40(%rcx),%xmm14

289
290	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
291	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
292	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000af7 <+2807>:	movaps 0x50(%r8),%xmm14
   0x0000000000000afc <+2812>:	mulps  0x50(%rcx),%xmm14
   0x0000000000000b01 <+2817>:	movss  %xmm3,(%rsp)
   0x0000000000000b06 <+2822>:	movaps 0x10(%r8),%xmm3
   0x0000000000000b0b <+2827>:	mulps  0x50(%rcx),%xmm3
   0x0000000000000b0f <+2831>:	movss  %xmm2,0x8(%rsp)
   0x0000000000000b24 <+2852>:	addps  %xmm15,%xmm3
   0x0000000000000b3b <+2875>:	addps  %xmm2,%xmm14
   0x0000000000000b7b <+2939>:	movaps 0x90(%r8),%xmm3
   0x0000000000000b83 <+2947>:	mulps  0x50(%rcx),%xmm3
   0x0000000000000b87 <+2951>:	addps  %xmm2,%xmm3
   0x0000000000000b8a <+2954>:	movaps 0xd0(%r8),%xmm2
   0x0000000000000b92 <+2962>:	mulps  0x50(%rcx),%xmm2
   0x0000000000000bb7 <+2999>:	addps  %xmm14,%xmm2

293
294	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
295	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
296	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000b28 <+2856>:	movaps 0x20(%r8),%xmm15
   0x0000000000000b32 <+2866>:	mulps  0x60(%rcx),%xmm15
   0x0000000000000b4b <+2891>:	addps  %xmm3,%xmm15
   0x0000000000000b4f <+2895>:	movaps 0x60(%r8),%xmm3
   0x0000000000000b54 <+2900>:	mulps  0x60(%rcx),%xmm3
   0x0000000000000b66 <+2918>:	addps  %xmm14,%xmm3
   0x0000000000000b6a <+2922>:	movaps 0xa0(%r8),%xmm14
   0x0000000000000b72 <+2930>:	mulps  0x60(%rcx),%xmm14
   0x0000000000000b96 <+2966>:	addps  %xmm3,%xmm14
   0x0000000000000bbb <+3003>:	movaps 0xe0(%r8),%xmm14
   0x0000000000000bc3 <+3011>:	mulps  0x60(%rcx),%xmm14
   0x0000000000000bc8 <+3016>:	addps  %xmm2,%xmm14

297
298	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
299	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
300	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000b2d <+2861>:	movaps 0x30(%r8),%xmm0
   0x0000000000000b37 <+2871>:	mulps  0x70(%rcx),%xmm0
   0x0000000000000b58 <+2904>:	addps  %xmm15,%xmm0
   0x0000000000000b5c <+2908>:	movaps 0x70(%r8),%xmm15
   0x0000000000000b61 <+2913>:	mulps  0x70(%rcx),%xmm15
   0x0000000000000b77 <+2935>:	addps  %xmm3,%xmm15
   0x0000000000000b9a <+2970>:	movaps 0xb0(%r8),%xmm3
   0x0000000000000ba2 <+2978>:	mulps  0x70(%rcx),%xmm3
   0x0000000000000ba6 <+2982>:	addps  %xmm14,%xmm3
   0x0000000000000bcc <+3020>:	movaps 0xf0(%r8),%xmm2
   0x0000000000000bd4 <+3028>:	mulps  0x70(%rcx),%xmm2
   0x0000000000000bd8 <+3032>:	addps  %xmm14,%xmm2

301	      }
302
303	      /* A(1,2)*B(2,2) = C(1,2). */
304	      for (i = 0; i < 4; i++)
305	      {
306	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
307	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
308	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000bdc <+3036>:	movaps 0x100(%r8),%xmm14
   0x0000000000000be4 <+3044>:	mulps  0x140(%rcx),%xmm14
   0x0000000000000bec <+3052>:	addps  %xmm0,%xmm14
   0x0000000000000c2a <+3114>:	movaps 0x140(%r8),%xmm14
   0x0000000000000c32 <+3122>:	mulps  0x140(%rcx),%xmm14
   0x0000000000000c3a <+3130>:	addps  %xmm15,%xmm14
   0x0000000000000c7a <+3194>:	movaps 0x180(%r8),%xmm14
   0x0000000000000c82 <+3202>:	mulps  0x140(%rcx),%xmm14
   0x0000000000000c8a <+3210>:	addps  %xmm3,%xmm14
   0x0000000000000cc8 <+3272>:	movaps 0x1c0(%r8),%xmm14
   0x0000000000000cd0 <+3280>:	mulps  0x140(%rcx),%xmm14
   0x0000000000000cd8 <+3288>:	addps  %xmm2,%xmm14

309
310	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
311	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
312	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000bf0 <+3056>:	movaps 0x110(%r8),%xmm0
   0x0000000000000bf8 <+3064>:	mulps  0x150(%rcx),%xmm0
   0x0000000000000bff <+3071>:	addps  %xmm14,%xmm0
   0x0000000000000c3e <+3134>:	movaps 0x150(%r8),%xmm15
   0x0000000000000c46 <+3142>:	mulps  0x150(%rcx),%xmm15
   0x0000000000000c4e <+3150>:	addps  %xmm14,%xmm15
   0x0000000000000c8e <+3214>:	movaps 0x190(%r8),%xmm3
   0x0000000000000c96 <+3222>:	mulps  0x150(%rcx),%xmm3
   0x0000000000000c9d <+3229>:	addps  %xmm14,%xmm3
   0x0000000000000cdc <+3292>:	movaps 0x1d0(%r8),%xmm2
   0x0000000000000ce4 <+3300>:	mulps  0x150(%rcx),%xmm2
   0x0000000000000ceb <+3307>:	addps  %xmm14,%xmm2

313
314	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
315	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
316	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000c03 <+3075>:	movaps 0x120(%r8),%xmm14
   0x0000000000000c0b <+3083>:	mulps  0x160(%rcx),%xmm14
   0x0000000000000c13 <+3091>:	addps  %xmm0,%xmm14
   0x0000000000000c52 <+3154>:	movaps 0x160(%r8),%xmm14
   0x0000000000000c5a <+3162>:	mulps  0x160(%rcx),%xmm14
   0x0000000000000c62 <+3170>:	addps  %xmm15,%xmm14
   0x0000000000000ca1 <+3233>:	movaps 0x1a0(%r8),%xmm14
   0x0000000000000ca9 <+3241>:	mulps  0x160(%rcx),%xmm14
   0x0000000000000cb1 <+3249>:	addps  %xmm3,%xmm14
   0x0000000000000cef <+3311>:	movaps 0x1e0(%r8),%xmm14
   0x0000000000000cf7 <+3319>:	mulps  0x160(%rcx),%xmm14
   0x0000000000000cff <+3327>:	addps  %xmm2,%xmm14

317
318	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
319	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
320	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000c17 <+3095>:	movaps 0x130(%r8),%xmm0
   0x0000000000000c1f <+3103>:	mulps  0x170(%rcx),%xmm0
   0x0000000000000c26 <+3110>:	addps  %xmm14,%xmm0
   0x0000000000000c66 <+3174>:	movaps 0x170(%r8),%xmm15
   0x0000000000000c6e <+3182>:	mulps  0x170(%rcx),%xmm15
   0x0000000000000c76 <+3190>:	addps  %xmm14,%xmm15
   0x0000000000000cb5 <+3253>:	movaps 0x1b0(%r8),%xmm3
   0x0000000000000cbd <+3261>:	mulps  0x170(%rcx),%xmm3
   0x0000000000000cc4 <+3268>:	addps  %xmm14,%xmm3
   0x0000000000000d03 <+3331>:	movaps 0x1f0(%r8),%xmm2
   0x0000000000000d0b <+3339>:	mulps  0x170(%rcx),%xmm2
   0x0000000000000d12 <+3346>:	addps  %xmm14,%xmm2

321	      }
322
323	      /* A(1,3)*B(3,2) = C(1,2). */
324	      for (i = 0; i < 4; i++)
325	      {
326	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
327	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
328	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d16 <+3350>:	movaps 0x200(%r8),%xmm14
   0x0000000000000d1e <+3358>:	mulps  0x240(%rcx),%xmm14
   0x0000000000000d26 <+3366>:	addps  %xmm0,%xmm14
   0x0000000000000d64 <+3428>:	movaps 0x240(%r8),%xmm14
   0x0000000000000d6c <+3436>:	mulps  0x240(%rcx),%xmm14
   0x0000000000000d74 <+3444>:	addps  %xmm15,%xmm14
   0x0000000000000db4 <+3508>:	movaps 0x280(%r8),%xmm14
   0x0000000000000dbc <+3516>:	mulps  0x240(%rcx),%xmm14
   0x0000000000000dc4 <+3524>:	addps  %xmm3,%xmm14
   0x0000000000000e02 <+3586>:	movaps 0x2c0(%r8),%xmm14
   0x0000000000000e0a <+3594>:	mulps  0x240(%rcx),%xmm14
   0x0000000000000e12 <+3602>:	addps  %xmm2,%xmm14

329
330	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
331	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
332	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d2a <+3370>:	movaps 0x210(%r8),%xmm0
   0x0000000000000d32 <+3378>:	mulps  0x250(%rcx),%xmm0
   0x0000000000000d39 <+3385>:	addps  %xmm14,%xmm0
   0x0000000000000d78 <+3448>:	movaps 0x250(%r8),%xmm15
   0x0000000000000d80 <+3456>:	mulps  0x250(%rcx),%xmm15
   0x0000000000000d88 <+3464>:	addps  %xmm14,%xmm15
   0x0000000000000dc8 <+3528>:	movaps 0x290(%r8),%xmm3
   0x0000000000000dd0 <+3536>:	mulps  0x250(%rcx),%xmm3
   0x0000000000000dd7 <+3543>:	addps  %xmm14,%xmm3
   0x0000000000000e16 <+3606>:	movaps 0x2d0(%r8),%xmm2
   0x0000000000000e1e <+3614>:	mulps  0x250(%rcx),%xmm2
   0x0000000000000e25 <+3621>:	addps  %xmm14,%xmm2

333
334	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
335	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
336	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d3d <+3389>:	movaps 0x220(%r8),%xmm14
   0x0000000000000d45 <+3397>:	mulps  0x260(%rcx),%xmm14
   0x0000000000000d4d <+3405>:	addps  %xmm0,%xmm14
   0x0000000000000d8c <+3468>:	movaps 0x260(%r8),%xmm14
   0x0000000000000d94 <+3476>:	mulps  0x260(%rcx),%xmm14
   0x0000000000000d9c <+3484>:	addps  %xmm15,%xmm14
   0x0000000000000ddb <+3547>:	movaps 0x2a0(%r8),%xmm14
   0x0000000000000de3 <+3555>:	mulps  0x260(%rcx),%xmm14
   0x0000000000000deb <+3563>:	addps  %xmm3,%xmm14
   0x0000000000000e29 <+3625>:	movaps 0x2e0(%r8),%xmm14
   0x0000000000000e31 <+3633>:	mulps  0x260(%rcx),%xmm14
   0x0000000000000e39 <+3641>:	addps  %xmm2,%xmm14

337
338	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
339	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
340	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d51 <+3409>:	movaps 0x230(%r8),%xmm0
   0x0000000000000d59 <+3417>:	mulps  0x270(%rcx),%xmm0
   0x0000000000000d60 <+3424>:	addps  %xmm14,%xmm0
   0x0000000000000da0 <+3488>:	movaps 0x270(%r8),%xmm15
   0x0000000000000da8 <+3496>:	mulps  0x270(%rcx),%xmm15
   0x0000000000000db0 <+3504>:	addps  %xmm14,%xmm15
   0x0000000000000def <+3567>:	movaps 0x2b0(%r8),%xmm3
   0x0000000000000df7 <+3575>:	mulps  0x270(%rcx),%xmm3
   0x0000000000000dfe <+3582>:	addps  %xmm14,%xmm3
   0x0000000000000e3d <+3645>:	movaps 0x2f0(%r8),%xmm2
   0x0000000000000e45 <+3653>:	mulps  0x270(%rcx),%xmm2
   0x0000000000000e4c <+3660>:	addps  %xmm14,%xmm2

341	      }
342
343	      /* A(1,4)*B(4,2) = C(1,2). */
344	      for (i = 0; i < 4; i++)
345	      {
346	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
347	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
348	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000e50 <+3664>:	movaps 0x300(%r8),%xmm14
   0x0000000000000e58 <+3672>:	mulps  0x340(%rcx),%xmm14
   0x0000000000000e60 <+3680>:	addps  %xmm0,%xmm14
   0x0000000000000e9e <+3742>:	movaps 0x340(%r8),%xmm14
   0x0000000000000ea6 <+3750>:	mulps  0x340(%rcx),%xmm14
   0x0000000000000eae <+3758>:	addps  %xmm15,%xmm14
   0x0000000000000eee <+3822>:	movaps 0x380(%r8),%xmm14
   0x0000000000000ef6 <+3830>:	mulps  0x340(%rcx),%xmm14
   0x0000000000000efe <+3838>:	addps  %xmm3,%xmm14
   0x0000000000000f3c <+3900>:	movaps 0x3c0(%r8),%xmm14
   0x0000000000000f44 <+3908>:	mulps  0x340(%rcx),%xmm14
   0x0000000000000f4c <+3916>:	addps  %xmm2,%xmm14

349
350	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
351	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
352	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000e64 <+3684>:	movaps 0x310(%r8),%xmm0
   0x0000000000000e6c <+3692>:	mulps  0x350(%rcx),%xmm0
   0x0000000000000e73 <+3699>:	addps  %xmm14,%xmm0
   0x0000000000000eb2 <+3762>:	movaps 0x350(%r8),%xmm15
   0x0000000000000eba <+3770>:	mulps  0x350(%rcx),%xmm15
   0x0000000000000ec2 <+3778>:	addps  %xmm14,%xmm15
   0x0000000000000f02 <+3842>:	movaps 0x390(%r8),%xmm3
   0x0000000000000f0a <+3850>:	mulps  0x350(%rcx),%xmm3
   0x0000000000000f11 <+3857>:	addps  %xmm14,%xmm3
   0x0000000000000f50 <+3920>:	movaps 0x3d0(%r8),%xmm2
   0x0000000000000f58 <+3928>:	mulps  0x350(%rcx),%xmm2
   0x0000000000000f5f <+3935>:	addps  %xmm14,%xmm2

353
354	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
355	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
356	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000e77 <+3703>:	movaps 0x320(%r8),%xmm14
   0x0000000000000e7f <+3711>:	mulps  0x360(%rcx),%xmm14
   0x0000000000000e87 <+3719>:	addps  %xmm0,%xmm14
   0x0000000000000ec6 <+3782>:	movaps 0x360(%r8),%xmm14
   0x0000000000000ece <+3790>:	mulps  0x360(%rcx),%xmm14
   0x0000000000000ed6 <+3798>:	addps  %xmm15,%xmm14
   0x0000000000000f15 <+3861>:	movaps 0x3a0(%r8),%xmm14
   0x0000000000000f1d <+3869>:	mulps  0x360(%rcx),%xmm14
   0x0000000000000f25 <+3877>:	addps  %xmm3,%xmm14
   0x0000000000000f63 <+3939>:	movaps 0x3e0(%r8),%xmm14
   0x0000000000000f6b <+3947>:	mulps  0x360(%rcx),%xmm14
   0x0000000000000f73 <+3955>:	addps  %xmm2,%xmm14

357
358	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
359	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
360	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000e8b <+3723>:	movaps 0x330(%r8),%xmm0
   0x0000000000000e93 <+3731>:	mulps  0x370(%rcx),%xmm0
   0x0000000000000e9a <+3738>:	addps  %xmm14,%xmm0
   0x0000000000000eda <+3802>:	movaps 0x370(%r8),%xmm15
   0x0000000000000ee2 <+3810>:	mulps  0x370(%rcx),%xmm15
   0x0000000000000eea <+3818>:	addps  %xmm14,%xmm15
   0x0000000000000f29 <+3881>:	movaps 0x3b0(%r8),%xmm3
   0x0000000000000f31 <+3889>:	mulps  0x370(%rcx),%xmm3
   0x0000000000000f38 <+3896>:	addps  %xmm14,%xmm3
   0x0000000000000f77 <+3959>:	movaps 0x3f0(%r8),%xmm2
   0x0000000000000f7f <+3967>:	mulps  0x370(%rcx),%xmm2
   0x0000000000000f86 <+3974>:	addps  %xmm14,%xmm2

361	      }
362
363	      /* Store C(1,2) block. */
364	      for (i = 0; i < 4; i++)
365	      {
366	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000000f8a <+3978>:	movaps 0x20(%rsp),%xmm14
   0x0000000000000f90 <+3984>:	mulps  %xmm14,%xmm0
   0x0000000000000fa4 <+4004>:	mulps  %xmm14,%xmm15
   0x0000000000000fb2 <+4018>:	mulps  %xmm14,%xmm3
   0x0000000000000fc5 <+4037>:	mulps  %xmm14,%xmm2

367	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_12]), C_row[i]);
   0x0000000000000f94 <+3988>:	addps  0x40(%r9),%xmm0
   0x0000000000000fa8 <+4008>:	addps  0x50(%r9),%xmm15
   0x0000000000000fb6 <+4022>:	addps  0x60(%r9),%xmm3
   0x0000000000000fc9 <+4041>:	addps  0x70(%r9),%xmm2

368	        _mm_store_ps(&C[i*4+C_OFFSET_12], C_row[i]);
   0x0000000000000f99 <+3993>:	movaps %xmm0,0x40(%r9)
   0x0000000000000f9e <+3998>:	movss  0x10(%rsp),%xmm0
   0x0000000000000fad <+4013>:	movaps %xmm15,0x50(%r9)
   0x0000000000000fbb <+4027>:	movaps %xmm3,0x60(%r9)
   0x0000000000000fc0 <+4032>:	movss  (%rsp),%xmm3
   0x0000000000000fce <+4046>:	movaps %xmm2,0x70(%r9)
   0x0000000000000fd3 <+4051>:	movss  0x8(%rsp),%xmm2

369	      }
370	    }
371
372	    /* Reset C(1,3) matrix accumulators */
373	    C_row[0] = _mm_setzero_ps();
374	    C_row[1] = _mm_setzero_ps();
375	    C_row[2] = _mm_setzero_ps();
376	    C_row[3] = _mm_setzero_ps();
377
378	    if (norm_product[0][2] &&
   0x0000000000000fd9 <+4057>:	comiss %xmm1,%xmm13
   0x0000000000000fdd <+4061>:	jb     0x1553 <stream_kernel+5459>

379	        norm_product[1][6] &&
   0x0000000000000fe3 <+4067>:	movss  0x170(%rsp),%xmm13
   0x0000000000000fed <+4077>:	comiss %xmm1,%xmm13
   0x0000000000000ff1 <+4081>:	jb     0x1553 <stream_kernel+5459>

380	        norm_product[2][10] &&
   0x0000000000000ff7 <+4087>:	movss  0x168(%rsp),%xmm13
   0x0000000000001001 <+4097>:	comiss %xmm1,%xmm13
   0x0000000000001005 <+4101>:	jb     0x1553 <stream_kernel+5459>

381	        norm_product[3][14])
   0x000000000000100b <+4107>:	movss  0x158(%rsp),%xmm13
   0x0000000000001015 <+4117>:	comiss %xmm1,%xmm13
   0x0000000000001019 <+4121>:	jb     0x1553 <stream_kernel+5459>

382	    {
383	      /* A(1,1)*B(1,3) = C(1,3). */
384	      for (i = 0; i < 4; i++)
385	      {
386	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
387	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
388	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000101f <+4127>:	movaps (%r8),%xmm15
   0x000000000000102d <+4141>:	mulps  0x80(%rcx),%xmm15
   0x000000000000104b <+4171>:	movaps 0x40(%r8),%xmm2
   0x0000000000001050 <+4176>:	mulps  0x80(%rcx),%xmm2
   0x0000000000001057 <+4183>:	movss  %xmm0,0x10(%rsp)
   0x000000000000107e <+4222>:	movaps 0x80(%r8),%xmm2
   0x0000000000001086 <+4230>:	mulps  0x80(%rcx),%xmm2
   0x0000000000001102 <+4354>:	movaps 0xc0(%r8),%xmm14
   0x000000000000110a <+4362>:	mulps  0x80(%rcx),%xmm14

389
390	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
391	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
392	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001023 <+4131>:	movaps 0x10(%r8),%xmm13
   0x0000000000001028 <+4136>:	movaps 0x50(%r8),%xmm14
   0x0000000000001035 <+4149>:	mulps  0x90(%rcx),%xmm13
   0x000000000000103d <+4157>:	mulps  0x90(%rcx),%xmm14
   0x0000000000001045 <+4165>:	movss  %xmm2,0x8(%rsp)
   0x0000000000001062 <+4194>:	addps  %xmm15,%xmm13
   0x000000000000107a <+4218>:	addps  %xmm2,%xmm14
   0x00000000000010c7 <+4295>:	movaps 0x90(%r8),%xmm13
   0x00000000000010cf <+4303>:	mulps  0x90(%rcx),%xmm13
   0x00000000000010d7 <+4311>:	addps  %xmm2,%xmm13
   0x00000000000010db <+4315>:	movaps 0xd0(%r8),%xmm2
   0x00000000000010e3 <+4323>:	mulps  0x90(%rcx),%xmm2
   0x0000000000001112 <+4370>:	addps  %xmm14,%xmm2

393
394	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
395	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
396	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001066 <+4198>:	movaps 0x20(%r8),%xmm15
   0x000000000000106b <+4203>:	mulps  0xa0(%rcx),%xmm15
   0x000000000000108d <+4237>:	addps  %xmm13,%xmm15
   0x0000000000001091 <+4241>:	movaps 0x60(%r8),%xmm13
   0x0000000000001096 <+4246>:	mulps  0xa0(%rcx),%xmm13
   0x00000000000010af <+4271>:	addps  %xmm14,%xmm13
   0x00000000000010b3 <+4275>:	movaps 0xa0(%r8),%xmm14
   0x00000000000010bb <+4283>:	mulps  0xa0(%rcx),%xmm14
   0x00000000000010ea <+4330>:	addps  %xmm13,%xmm14
   0x0000000000001116 <+4374>:	movaps 0xe0(%r8),%xmm14
   0x000000000000111e <+4382>:	mulps  0xa0(%rcx),%xmm14
   0x0000000000001126 <+4390>:	addps  %xmm2,%xmm14

397
398	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
399	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
400	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000105d <+4189>:	movaps 0x30(%r8),%xmm0
   0x0000000000001073 <+4211>:	mulps  0xb0(%rcx),%xmm0
   0x000000000000109e <+4254>:	addps  %xmm15,%xmm0
   0x00000000000010a2 <+4258>:	movaps 0x70(%r8),%xmm15
   0x00000000000010a7 <+4263>:	mulps  0xb0(%rcx),%xmm15
   0x00000000000010c3 <+4291>:	addps  %xmm13,%xmm15
   0x00000000000010ee <+4334>:	movaps 0xb0(%r8),%xmm13
   0x00000000000010f6 <+4342>:	mulps  0xb0(%rcx),%xmm13
   0x00000000000010fe <+4350>:	addps  %xmm14,%xmm13
   0x000000000000112a <+4394>:	movaps 0xf0(%r8),%xmm2
   0x0000000000001132 <+4402>:	mulps  0xb0(%rcx),%xmm2
   0x0000000000001139 <+4409>:	addps  %xmm14,%xmm2

401	      }
402
403	      /* A(1,2)*B(2,3) = C(1,3). */
404	      for (i = 0; i < 4; i++)
405	      {
406	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
407	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
408	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000113d <+4413>:	movaps 0x100(%r8),%xmm14
   0x0000000000001145 <+4421>:	mulps  0x180(%rcx),%xmm14
   0x000000000000114d <+4429>:	addps  %xmm0,%xmm14
   0x000000000000118b <+4491>:	movaps 0x140(%r8),%xmm14
   0x0000000000001193 <+4499>:	mulps  0x180(%rcx),%xmm14
   0x000000000000119b <+4507>:	addps  %xmm15,%xmm14
   0x00000000000011db <+4571>:	movaps 0x180(%r8),%xmm14
   0x00000000000011e3 <+4579>:	mulps  0x180(%rcx),%xmm14
   0x00000000000011eb <+4587>:	addps  %xmm13,%xmm14
   0x000000000000122b <+4651>:	movaps 0x1c0(%r8),%xmm14
   0x0000000000001233 <+4659>:	mulps  0x180(%rcx),%xmm14
   0x000000000000123b <+4667>:	addps  %xmm2,%xmm14

409
410	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
411	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
412	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001151 <+4433>:	movaps 0x110(%r8),%xmm0
   0x0000000000001159 <+4441>:	mulps  0x190(%rcx),%xmm0
   0x0000000000001160 <+4448>:	addps  %xmm14,%xmm0
   0x000000000000119f <+4511>:	movaps 0x150(%r8),%xmm15
   0x00000000000011a7 <+4519>:	mulps  0x190(%rcx),%xmm15
   0x00000000000011af <+4527>:	addps  %xmm14,%xmm15
   0x00000000000011ef <+4591>:	movaps 0x190(%r8),%xmm13
   0x00000000000011f7 <+4599>:	mulps  0x190(%rcx),%xmm13
   0x00000000000011ff <+4607>:	addps  %xmm14,%xmm13
   0x000000000000123f <+4671>:	movaps 0x1d0(%r8),%xmm2
   0x0000000000001247 <+4679>:	mulps  0x190(%rcx),%xmm2
   0x000000000000124e <+4686>:	addps  %xmm14,%xmm2

413
414	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
415	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
416	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001164 <+4452>:	movaps 0x120(%r8),%xmm14
   0x000000000000116c <+4460>:	mulps  0x1a0(%rcx),%xmm14
   0x0000000000001174 <+4468>:	addps  %xmm0,%xmm14
   0x00000000000011b3 <+4531>:	movaps 0x160(%r8),%xmm14
   0x00000000000011bb <+4539>:	mulps  0x1a0(%rcx),%xmm14
   0x00000000000011c3 <+4547>:	addps  %xmm15,%xmm14
   0x0000000000001203 <+4611>:	movaps 0x1a0(%r8),%xmm14
   0x000000000000120b <+4619>:	mulps  0x1a0(%rcx),%xmm14
   0x0000000000001213 <+4627>:	addps  %xmm13,%xmm14
   0x0000000000001252 <+4690>:	movaps 0x1e0(%r8),%xmm14
   0x000000000000125a <+4698>:	mulps  0x1a0(%rcx),%xmm14
   0x0000000000001262 <+4706>:	addps  %xmm2,%xmm14

417
418	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
419	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
420	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001178 <+4472>:	movaps 0x130(%r8),%xmm0
   0x0000000000001180 <+4480>:	mulps  0x1b0(%rcx),%xmm0
   0x0000000000001187 <+4487>:	addps  %xmm14,%xmm0
   0x00000000000011c7 <+4551>:	movaps 0x170(%r8),%xmm15
   0x00000000000011cf <+4559>:	mulps  0x1b0(%rcx),%xmm15
   0x00000000000011d7 <+4567>:	addps  %xmm14,%xmm15
   0x0000000000001217 <+4631>:	movaps 0x1b0(%r8),%xmm13
   0x000000000000121f <+4639>:	mulps  0x1b0(%rcx),%xmm13
   0x0000000000001227 <+4647>:	addps  %xmm14,%xmm13
   0x0000000000001266 <+4710>:	movaps 0x1f0(%r8),%xmm2
   0x000000000000126e <+4718>:	mulps  0x1b0(%rcx),%xmm2
   0x0000000000001275 <+4725>:	addps  %xmm14,%xmm2

421	      }
422
423	      /* A(1,3)*B(3,3) = C(1,3). */
424	      for (i = 0; i < 4; i++)
425	      {
426	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
427	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
428	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001279 <+4729>:	movaps 0x200(%r8),%xmm14
   0x0000000000001281 <+4737>:	mulps  0x280(%rcx),%xmm14
   0x0000000000001289 <+4745>:	addps  %xmm0,%xmm14
   0x00000000000012c7 <+4807>:	movaps 0x240(%r8),%xmm14
   0x00000000000012cf <+4815>:	mulps  0x280(%rcx),%xmm14
   0x00000000000012d7 <+4823>:	addps  %xmm15,%xmm14
   0x0000000000001317 <+4887>:	movaps 0x280(%r8),%xmm14
   0x000000000000131f <+4895>:	mulps  0x280(%rcx),%xmm14
   0x0000000000001327 <+4903>:	addps  %xmm13,%xmm14
   0x0000000000001367 <+4967>:	movaps 0x2c0(%r8),%xmm14
   0x000000000000136f <+4975>:	mulps  0x280(%rcx),%xmm14
   0x0000000000001377 <+4983>:	addps  %xmm2,%xmm14

429
430	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
431	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
432	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000128d <+4749>:	movaps 0x210(%r8),%xmm0
   0x0000000000001295 <+4757>:	mulps  0x290(%rcx),%xmm0
   0x000000000000129c <+4764>:	addps  %xmm14,%xmm0
   0x00000000000012db <+4827>:	movaps 0x250(%r8),%xmm15
   0x00000000000012e3 <+4835>:	mulps  0x290(%rcx),%xmm15
   0x00000000000012eb <+4843>:	addps  %xmm14,%xmm15
   0x000000000000132b <+4907>:	movaps 0x290(%r8),%xmm13
   0x0000000000001333 <+4915>:	mulps  0x290(%rcx),%xmm13
   0x000000000000133b <+4923>:	addps  %xmm14,%xmm13
   0x000000000000137b <+4987>:	movaps 0x2d0(%r8),%xmm2
   0x0000000000001383 <+4995>:	mulps  0x290(%rcx),%xmm2
   0x000000000000138a <+5002>:	addps  %xmm14,%xmm2

433
434	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
435	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
436	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000012a0 <+4768>:	movaps 0x220(%r8),%xmm14
   0x00000000000012a8 <+4776>:	mulps  0x2a0(%rcx),%xmm14
   0x00000000000012b0 <+4784>:	addps  %xmm0,%xmm14
   0x00000000000012ef <+4847>:	movaps 0x260(%r8),%xmm14
   0x00000000000012f7 <+4855>:	mulps  0x2a0(%rcx),%xmm14
   0x00000000000012ff <+4863>:	addps  %xmm15,%xmm14
   0x000000000000133f <+4927>:	movaps 0x2a0(%r8),%xmm14
   0x0000000000001347 <+4935>:	mulps  0x2a0(%rcx),%xmm14
   0x000000000000134f <+4943>:	addps  %xmm13,%xmm14
   0x000000000000138e <+5006>:	movaps 0x2e0(%r8),%xmm14
   0x0000000000001396 <+5014>:	mulps  0x2a0(%rcx),%xmm14
   0x000000000000139e <+5022>:	addps  %xmm2,%xmm14

437
438	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
439	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
440	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000012b4 <+4788>:	movaps 0x230(%r8),%xmm0
   0x00000000000012bc <+4796>:	mulps  0x2b0(%rcx),%xmm0
   0x00000000000012c3 <+4803>:	addps  %xmm14,%xmm0
   0x0000000000001303 <+4867>:	movaps 0x270(%r8),%xmm15
   0x000000000000130b <+4875>:	mulps  0x2b0(%rcx),%xmm15
   0x0000000000001313 <+4883>:	addps  %xmm14,%xmm15
   0x0000000000001353 <+4947>:	movaps 0x2b0(%r8),%xmm13
   0x000000000000135b <+4955>:	mulps  0x2b0(%rcx),%xmm13
   0x0000000000001363 <+4963>:	addps  %xmm14,%xmm13
   0x00000000000013a2 <+5026>:	movaps 0x2f0(%r8),%xmm2
   0x00000000000013aa <+5034>:	mulps  0x2b0(%rcx),%xmm2
   0x00000000000013b1 <+5041>:	addps  %xmm14,%xmm2

441	      }
442
443	      /* A(1,4)*B(4,3) = C(1,3). */
444	      for (i = 0; i < 4; i++)
445	      {
446	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
447	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
448	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000013b5 <+5045>:	movaps 0x300(%r8),%xmm14
   0x00000000000013bd <+5053>:	mulps  0x380(%rcx),%xmm14
   0x00000000000013c5 <+5061>:	addps  %xmm0,%xmm14
   0x0000000000001403 <+5123>:	movaps 0x340(%r8),%xmm14
   0x000000000000140b <+5131>:	mulps  0x380(%rcx),%xmm14
   0x0000000000001413 <+5139>:	addps  %xmm15,%xmm14
   0x0000000000001453 <+5203>:	movaps 0x380(%r8),%xmm14
   0x000000000000145b <+5211>:	mulps  0x380(%rcx),%xmm14
   0x0000000000001463 <+5219>:	addps  %xmm13,%xmm14
   0x00000000000014a3 <+5283>:	movaps 0x3c0(%r8),%xmm14
   0x00000000000014ab <+5291>:	mulps  0x380(%rcx),%xmm14
   0x00000000000014b3 <+5299>:	addps  %xmm2,%xmm14

449
450	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
451	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
452	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000013c9 <+5065>:	movaps 0x310(%r8),%xmm0
   0x00000000000013d1 <+5073>:	mulps  0x390(%rcx),%xmm0
   0x00000000000013d8 <+5080>:	addps  %xmm14,%xmm0
   0x0000000000001417 <+5143>:	movaps 0x350(%r8),%xmm15
   0x000000000000141f <+5151>:	mulps  0x390(%rcx),%xmm15
   0x0000000000001427 <+5159>:	addps  %xmm14,%xmm15
   0x0000000000001467 <+5223>:	movaps 0x390(%r8),%xmm13
   0x000000000000146f <+5231>:	mulps  0x390(%rcx),%xmm13
   0x0000000000001477 <+5239>:	addps  %xmm14,%xmm13
   0x00000000000014b7 <+5303>:	movaps 0x3d0(%r8),%xmm2
   0x00000000000014bf <+5311>:	mulps  0x390(%rcx),%xmm2
   0x00000000000014c6 <+5318>:	addps  %xmm14,%xmm2

453
454	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
455	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
456	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000013dc <+5084>:	movaps 0x320(%r8),%xmm14
   0x00000000000013e4 <+5092>:	mulps  0x3a0(%rcx),%xmm14
   0x00000000000013ec <+5100>:	addps  %xmm0,%xmm14
   0x000000000000142b <+5163>:	movaps 0x360(%r8),%xmm14
   0x0000000000001433 <+5171>:	mulps  0x3a0(%rcx),%xmm14
   0x000000000000143b <+5179>:	addps  %xmm15,%xmm14
   0x000000000000147b <+5243>:	movaps 0x3a0(%r8),%xmm14
   0x0000000000001483 <+5251>:	mulps  0x3a0(%rcx),%xmm14
   0x000000000000148b <+5259>:	addps  %xmm13,%xmm14
   0x00000000000014ca <+5322>:	movaps 0x3e0(%r8),%xmm14
   0x00000000000014d2 <+5330>:	mulps  0x3a0(%rcx),%xmm14
   0x00000000000014da <+5338>:	addps  %xmm2,%xmm14

457
458	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
459	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
460	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000013f0 <+5104>:	movaps 0x330(%r8),%xmm0
   0x00000000000013f8 <+5112>:	mulps  0x3b0(%rcx),%xmm0
   0x00000000000013ff <+5119>:	addps  %xmm14,%xmm0
   0x000000000000143f <+5183>:	movaps 0x370(%r8),%xmm15
   0x0000000000001447 <+5191>:	mulps  0x3b0(%rcx),%xmm15
   0x000000000000144f <+5199>:	addps  %xmm14,%xmm15
   0x000000000000148f <+5263>:	movaps 0x3b0(%r8),%xmm13
   0x0000000000001497 <+5271>:	mulps  0x3b0(%rcx),%xmm13
   0x000000000000149f <+5279>:	addps  %xmm14,%xmm13
   0x00000000000014de <+5342>:	movaps 0x3f0(%r8),%xmm2
   0x00000000000014e6 <+5350>:	mulps  0x3b0(%rcx),%xmm2
   0x00000000000014ed <+5357>:	addps  %xmm14,%xmm2

461	      }
462
463	      /* Store C(1,3) block. */
464	      for (i = 0; i < 4; i++)
465	      {
466	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000014f1 <+5361>:	movaps 0x20(%rsp),%xmm14
   0x00000000000014f7 <+5367>:	mulps  %xmm14,%xmm0
   0x0000000000001511 <+5393>:	mulps  %xmm14,%xmm15
   0x0000000000001525 <+5413>:	mulps  %xmm14,%xmm13
   0x0000000000001539 <+5433>:	mulps  %xmm14,%xmm2

467	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_13]), C_row[i]);
   0x00000000000014fb <+5371>:	addps  0x80(%r9),%xmm0
   0x0000000000001515 <+5397>:	addps  0x90(%r9),%xmm15
   0x0000000000001529 <+5417>:	addps  0xa0(%r9),%xmm13
   0x000000000000153d <+5437>:	addps  0xb0(%r9),%xmm2

468	        _mm_store_ps(&C[i*4+C_OFFSET_13], C_row[i]);
   0x0000000000001503 <+5379>:	movaps %xmm0,0x80(%r9)
   0x000000000000150b <+5387>:	movss  0x10(%rsp),%xmm0
   0x000000000000151d <+5405>:	movaps %xmm15,0x90(%r9)
   0x0000000000001531 <+5425>:	movaps %xmm13,0xa0(%r9)
   0x0000000000001545 <+5445>:	movaps %xmm2,0xb0(%r9)
   0x000000000000154d <+5453>:	movss  0x8(%rsp),%xmm2

469	      }
470	    }
471
472	    /* Reset C(1,4) matrix accumulators */
473	    C_row[0] = _mm_setzero_ps();
474	    C_row[1] = _mm_setzero_ps();
475	    C_row[2] = _mm_setzero_ps();
476	    C_row[3] = _mm_setzero_ps();
477
478	    if (norm_product[0][3] &&
   0x0000000000001553 <+5459>:	comiss %xmm1,%xmm12
   0x0000000000001557 <+5463>:	jb     0x1acb <stream_kernel+6859>

479	        norm_product[1][7] &&
   0x000000000000155d <+5469>:	movss  0x150(%rsp),%xmm12
   0x0000000000001567 <+5479>:	comiss %xmm1,%xmm12
   0x000000000000156b <+5483>:	jb     0x1acb <stream_kernel+6859>

480	        norm_product[2][11] &&
   0x0000000000001571 <+5489>:	movss  0x140(%rsp),%xmm12
   0x000000000000157b <+5499>:	comiss %xmm1,%xmm12
   0x000000000000157f <+5503>:	jb     0x1acb <stream_kernel+6859>

481	        norm_product[3][15])
   0x0000000000001585 <+5509>:	movss  0x138(%rsp),%xmm12
   0x000000000000158f <+5519>:	comiss %xmm1,%xmm12
   0x0000000000001593 <+5523>:	jb     0x1acb <stream_kernel+6859>

482	    {
483	      /* A(1,1)*B(1,4) = C(1,4). */
484	      for (i = 0; i < 4; i++)
485	      {
486	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
487	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
488	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001599 <+5529>:	movaps (%r8),%xmm15
   0x00000000000015a2 <+5538>:	movaps 0x40(%r8),%xmm12
   0x00000000000015ac <+5548>:	mulps  0xc0(%rcx),%xmm15
   0x00000000000015bc <+5564>:	mulps  0xc0(%rcx),%xmm12
   0x00000000000015f3 <+5619>:	movaps 0x80(%r8),%xmm12
   0x00000000000015fb <+5627>:	mulps  0xc0(%rcx),%xmm12
   0x0000000000001679 <+5753>:	movaps 0xc0(%r8),%xmm14
   0x0000000000001681 <+5761>:	mulps  0xc0(%rcx),%xmm14

489
490	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
491	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
492	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000159d <+5533>:	movaps 0x10(%r8),%xmm13
   0x00000000000015a7 <+5543>:	movaps 0x50(%r8),%xmm14
   0x00000000000015b4 <+5556>:	mulps  0xd0(%rcx),%xmm13
   0x00000000000015c4 <+5572>:	mulps  0xd0(%rcx),%xmm14
   0x00000000000015cc <+5580>:	movss  %xmm0,0x10(%rsp)
   0x00000000000015d7 <+5591>:	addps  %xmm15,%xmm13
   0x00000000000015ef <+5615>:	addps  %xmm12,%xmm14
   0x000000000000163d <+5693>:	movaps 0x90(%r8),%xmm13
   0x0000000000001645 <+5701>:	mulps  0xd0(%rcx),%xmm13
   0x000000000000164d <+5709>:	addps  %xmm12,%xmm13
   0x0000000000001651 <+5713>:	movaps 0xd0(%r8),%xmm12
   0x0000000000001659 <+5721>:	mulps  0xd0(%rcx),%xmm12
   0x0000000000001689 <+5769>:	addps  %xmm14,%xmm12

493
494	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
495	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
496	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000015db <+5595>:	movaps 0x20(%r8),%xmm15
   0x00000000000015e0 <+5600>:	mulps  0xe0(%rcx),%xmm15
   0x0000000000001603 <+5635>:	addps  %xmm13,%xmm15
   0x0000000000001607 <+5639>:	movaps 0x60(%r8),%xmm13
   0x000000000000160c <+5644>:	mulps  0xe0(%rcx),%xmm13
   0x0000000000001625 <+5669>:	addps  %xmm14,%xmm13
   0x0000000000001629 <+5673>:	movaps 0xa0(%r8),%xmm14
   0x0000000000001631 <+5681>:	mulps  0xe0(%rcx),%xmm14
   0x0000000000001661 <+5729>:	addps  %xmm13,%xmm14
   0x000000000000168d <+5773>:	movaps 0xe0(%r8),%xmm14
   0x0000000000001695 <+5781>:	mulps  0xe0(%rcx),%xmm14
   0x000000000000169d <+5789>:	addps  %xmm12,%xmm14

497
498	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
499	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
500	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000015d2 <+5586>:	movaps 0x30(%r8),%xmm0
   0x00000000000015e8 <+5608>:	mulps  0xf0(%rcx),%xmm0
   0x0000000000001614 <+5652>:	addps  %xmm15,%xmm0
   0x0000000000001618 <+5656>:	movaps 0x70(%r8),%xmm15
   0x000000000000161d <+5661>:	mulps  0xf0(%rcx),%xmm15
   0x0000000000001639 <+5689>:	addps  %xmm13,%xmm15
   0x0000000000001665 <+5733>:	movaps 0xb0(%r8),%xmm13
   0x000000000000166d <+5741>:	mulps  0xf0(%rcx),%xmm13
   0x0000000000001675 <+5749>:	addps  %xmm14,%xmm13
   0x00000000000016a1 <+5793>:	movaps 0xf0(%r8),%xmm12
   0x00000000000016a9 <+5801>:	mulps  0xf0(%rcx),%xmm12
   0x00000000000016b1 <+5809>:	addps  %xmm14,%xmm12

501	      }
502
503	      /* A(1,2)*B(2,4) = C(1,4). */
504	      for (i = 0; i < 4; i++)
505	      {
506	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
507	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
508	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000016b5 <+5813>:	movaps 0x100(%r8),%xmm14
   0x00000000000016bd <+5821>:	mulps  0x1c0(%rcx),%xmm14
   0x00000000000016c5 <+5829>:	addps  %xmm0,%xmm14
   0x0000000000001703 <+5891>:	movaps 0x140(%r8),%xmm14
   0x000000000000170b <+5899>:	mulps  0x1c0(%rcx),%xmm14
   0x0000000000001713 <+5907>:	addps  %xmm15,%xmm14
   0x0000000000001753 <+5971>:	movaps 0x180(%r8),%xmm14
   0x000000000000175b <+5979>:	mulps  0x1c0(%rcx),%xmm14
   0x0000000000001763 <+5987>:	addps  %xmm13,%xmm14
   0x00000000000017a3 <+6051>:	movaps 0x1c0(%r8),%xmm14
   0x00000000000017ab <+6059>:	mulps  0x1c0(%rcx),%xmm14
   0x00000000000017b3 <+6067>:	addps  %xmm12,%xmm14

509
510	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
511	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
512	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000016c9 <+5833>:	movaps 0x110(%r8),%xmm0
   0x00000000000016d1 <+5841>:	mulps  0x1d0(%rcx),%xmm0
   0x00000000000016d8 <+5848>:	addps  %xmm14,%xmm0
   0x0000000000001717 <+5911>:	movaps 0x150(%r8),%xmm15
   0x000000000000171f <+5919>:	mulps  0x1d0(%rcx),%xmm15
   0x0000000000001727 <+5927>:	addps  %xmm14,%xmm15
   0x0000000000001767 <+5991>:	movaps 0x190(%r8),%xmm13
   0x000000000000176f <+5999>:	mulps  0x1d0(%rcx),%xmm13
   0x0000000000001777 <+6007>:	addps  %xmm14,%xmm13
   0x00000000000017b7 <+6071>:	movaps 0x1d0(%r8),%xmm12
   0x00000000000017bf <+6079>:	mulps  0x1d0(%rcx),%xmm12
   0x00000000000017c7 <+6087>:	addps  %xmm14,%xmm12

513
514	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
515	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
516	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000016dc <+5852>:	movaps 0x120(%r8),%xmm14
   0x00000000000016e4 <+5860>:	mulps  0x1e0(%rcx),%xmm14
   0x00000000000016ec <+5868>:	addps  %xmm0,%xmm14
   0x000000000000172b <+5931>:	movaps 0x160(%r8),%xmm14
   0x0000000000001733 <+5939>:	mulps  0x1e0(%rcx),%xmm14
   0x000000000000173b <+5947>:	addps  %xmm15,%xmm14
   0x000000000000177b <+6011>:	movaps 0x1a0(%r8),%xmm14
   0x0000000000001783 <+6019>:	mulps  0x1e0(%rcx),%xmm14
   0x000000000000178b <+6027>:	addps  %xmm13,%xmm14
   0x00000000000017cb <+6091>:	movaps 0x1e0(%r8),%xmm14
   0x00000000000017d3 <+6099>:	mulps  0x1e0(%rcx),%xmm14
   0x00000000000017db <+6107>:	addps  %xmm12,%xmm14

517
518	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
519	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
520	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000016f0 <+5872>:	movaps 0x130(%r8),%xmm0
   0x00000000000016f8 <+5880>:	mulps  0x1f0(%rcx),%xmm0
   0x00000000000016ff <+5887>:	addps  %xmm14,%xmm0
   0x000000000000173f <+5951>:	movaps 0x170(%r8),%xmm15
   0x0000000000001747 <+5959>:	mulps  0x1f0(%rcx),%xmm15
   0x000000000000174f <+5967>:	addps  %xmm14,%xmm15
   0x000000000000178f <+6031>:	movaps 0x1b0(%r8),%xmm13
   0x0000000000001797 <+6039>:	mulps  0x1f0(%rcx),%xmm13
   0x000000000000179f <+6047>:	addps  %xmm14,%xmm13
   0x00000000000017df <+6111>:	movaps 0x1f0(%r8),%xmm12
   0x00000000000017e7 <+6119>:	mulps  0x1f0(%rcx),%xmm12
   0x00000000000017ef <+6127>:	addps  %xmm14,%xmm12

521	      }
522
523	      /* A(1,3)*B(3,4) = C(1,4). */
524	      for (i = 0; i < 4; i++)
525	      {
526	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
527	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
528	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000017f3 <+6131>:	movaps 0x200(%r8),%xmm14
   0x00000000000017fb <+6139>:	mulps  0x2c0(%rcx),%xmm14
   0x0000000000001803 <+6147>:	addps  %xmm0,%xmm14
   0x0000000000001841 <+6209>:	movaps 0x240(%r8),%xmm14
   0x0000000000001849 <+6217>:	mulps  0x2c0(%rcx),%xmm14
   0x0000000000001851 <+6225>:	addps  %xmm15,%xmm14
   0x0000000000001891 <+6289>:	movaps 0x280(%r8),%xmm14
   0x0000000000001899 <+6297>:	mulps  0x2c0(%rcx),%xmm14
   0x00000000000018a1 <+6305>:	addps  %xmm13,%xmm14
   0x00000000000018e1 <+6369>:	movaps 0x2c0(%r8),%xmm14
   0x00000000000018e9 <+6377>:	mulps  0x2c0(%rcx),%xmm14
   0x00000000000018f1 <+6385>:	addps  %xmm12,%xmm14

529
530	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
531	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
532	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001807 <+6151>:	movaps 0x210(%r8),%xmm0
   0x000000000000180f <+6159>:	mulps  0x2d0(%rcx),%xmm0
   0x0000000000001816 <+6166>:	addps  %xmm14,%xmm0
   0x0000000000001855 <+6229>:	movaps 0x250(%r8),%xmm15
   0x000000000000185d <+6237>:	mulps  0x2d0(%rcx),%xmm15
   0x0000000000001865 <+6245>:	addps  %xmm14,%xmm15
   0x00000000000018a5 <+6309>:	movaps 0x290(%r8),%xmm13
   0x00000000000018ad <+6317>:	mulps  0x2d0(%rcx),%xmm13
   0x00000000000018b5 <+6325>:	addps  %xmm14,%xmm13
   0x00000000000018f5 <+6389>:	movaps 0x2d0(%r8),%xmm12
   0x00000000000018fd <+6397>:	mulps  0x2d0(%rcx),%xmm12
   0x0000000000001905 <+6405>:	addps  %xmm14,%xmm12

533
534	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
535	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
536	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000181a <+6170>:	movaps 0x220(%r8),%xmm14
   0x0000000000001822 <+6178>:	mulps  0x2e0(%rcx),%xmm14
   0x000000000000182a <+6186>:	addps  %xmm0,%xmm14
   0x0000000000001869 <+6249>:	movaps 0x260(%r8),%xmm14
   0x0000000000001871 <+6257>:	mulps  0x2e0(%rcx),%xmm14
   0x0000000000001879 <+6265>:	addps  %xmm15,%xmm14
   0x00000000000018b9 <+6329>:	movaps 0x2a0(%r8),%xmm14
   0x00000000000018c1 <+6337>:	mulps  0x2e0(%rcx),%xmm14
   0x00000000000018c9 <+6345>:	addps  %xmm13,%xmm14
   0x0000000000001909 <+6409>:	movaps 0x2e0(%r8),%xmm14
   0x0000000000001911 <+6417>:	mulps  0x2e0(%rcx),%xmm14
   0x0000000000001919 <+6425>:	addps  %xmm12,%xmm14

537
538	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
539	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
540	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000182e <+6190>:	movaps 0x230(%r8),%xmm0
   0x0000000000001836 <+6198>:	mulps  0x2f0(%rcx),%xmm0
   0x000000000000183d <+6205>:	addps  %xmm14,%xmm0
   0x000000000000187d <+6269>:	movaps 0x270(%r8),%xmm15
   0x0000000000001885 <+6277>:	mulps  0x2f0(%rcx),%xmm15
   0x000000000000188d <+6285>:	addps  %xmm14,%xmm15
   0x00000000000018cd <+6349>:	movaps 0x2b0(%r8),%xmm13
   0x00000000000018d5 <+6357>:	mulps  0x2f0(%rcx),%xmm13
   0x00000000000018dd <+6365>:	addps  %xmm14,%xmm13
   0x000000000000191d <+6429>:	movaps 0x2f0(%r8),%xmm12
   0x0000000000001925 <+6437>:	mulps  0x2f0(%rcx),%xmm12
   0x000000000000192d <+6445>:	addps  %xmm14,%xmm12

541	      }
542
543	      /* A(1,4)*B(4,4) = C(1,4). */
544	      for (i = 0; i < 4; i++)
545	      {
546	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
547	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
548	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001931 <+6449>:	movaps 0x300(%r8),%xmm14
   0x0000000000001939 <+6457>:	mulps  0x3c0(%rcx),%xmm14
   0x0000000000001941 <+6465>:	addps  %xmm0,%xmm14
   0x000000000000197f <+6527>:	movaps 0x340(%r8),%xmm14
   0x0000000000001987 <+6535>:	mulps  0x3c0(%rcx),%xmm14
   0x000000000000198f <+6543>:	addps  %xmm15,%xmm14
   0x00000000000019cf <+6607>:	movaps 0x380(%r8),%xmm14
   0x00000000000019d7 <+6615>:	mulps  0x3c0(%rcx),%xmm14
   0x00000000000019df <+6623>:	addps  %xmm13,%xmm14
   0x0000000000001a1f <+6687>:	movaps 0x3c0(%r8),%xmm14
   0x0000000000001a27 <+6695>:	mulps  0x3c0(%rcx),%xmm14
   0x0000000000001a2f <+6703>:	addps  %xmm12,%xmm14

549
550	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
551	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
552	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001945 <+6469>:	movaps 0x310(%r8),%xmm0
   0x000000000000194d <+6477>:	mulps  0x3d0(%rcx),%xmm0
   0x0000000000001954 <+6484>:	addps  %xmm14,%xmm0
   0x0000000000001993 <+6547>:	movaps 0x350(%r8),%xmm15
   0x000000000000199b <+6555>:	mulps  0x3d0(%rcx),%xmm15
   0x00000000000019a3 <+6563>:	addps  %xmm14,%xmm15
   0x00000000000019e3 <+6627>:	movaps 0x390(%r8),%xmm13
   0x00000000000019eb <+6635>:	mulps  0x3d0(%rcx),%xmm13
   0x00000000000019f3 <+6643>:	addps  %xmm14,%xmm13
   0x0000000000001a33 <+6707>:	movaps 0x3d0(%r8),%xmm12
   0x0000000000001a3b <+6715>:	mulps  0x3d0(%rcx),%xmm12
   0x0000000000001a43 <+6723>:	addps  %xmm14,%xmm12

553
554	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
555	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
556	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001958 <+6488>:	movaps 0x320(%r8),%xmm14
   0x0000000000001960 <+6496>:	mulps  0x3e0(%rcx),%xmm14
   0x0000000000001968 <+6504>:	addps  %xmm0,%xmm14
   0x00000000000019a7 <+6567>:	movaps 0x360(%r8),%xmm14
   0x00000000000019af <+6575>:	mulps  0x3e0(%rcx),%xmm14
   0x00000000000019b7 <+6583>:	addps  %xmm15,%xmm14
   0x00000000000019f7 <+6647>:	movaps 0x3a0(%r8),%xmm14
   0x00000000000019ff <+6655>:	mulps  0x3e0(%rcx),%xmm14
   0x0000000000001a07 <+6663>:	addps  %xmm13,%xmm14
   0x0000000000001a47 <+6727>:	movaps 0x3e0(%r8),%xmm14
   0x0000000000001a4f <+6735>:	mulps  0x3e0(%rcx),%xmm14
   0x0000000000001a57 <+6743>:	addps  %xmm12,%xmm14

557
558	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
559	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
560	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000196c <+6508>:	movaps 0x330(%r8),%xmm0
   0x0000000000001974 <+6516>:	mulps  0x3f0(%rcx),%xmm0
   0x000000000000197b <+6523>:	addps  %xmm14,%xmm0
   0x00000000000019bb <+6587>:	movaps 0x370(%r8),%xmm15
   0x00000000000019c3 <+6595>:	mulps  0x3f0(%rcx),%xmm15
   0x00000000000019cb <+6603>:	addps  %xmm14,%xmm15
   0x0000000000001a0b <+6667>:	movaps 0x3b0(%r8),%xmm13
   0x0000000000001a13 <+6675>:	mulps  0x3f0(%rcx),%xmm13
   0x0000000000001a1b <+6683>:	addps  %xmm14,%xmm13
   0x0000000000001a5b <+6747>:	movaps 0x3f0(%r8),%xmm12
   0x0000000000001a63 <+6755>:	mulps  0x3f0(%rcx),%xmm12
   0x0000000000001a6b <+6763>:	addps  %xmm14,%xmm12

561	      }
562
563	      /* Store C(1,4) block. */
564	      for (i = 0; i < 4; i++)
565	      {
566	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000001a6f <+6767>:	movaps 0x20(%rsp),%xmm14
   0x0000000000001a75 <+6773>:	mulps  %xmm14,%xmm0
   0x0000000000001a8f <+6799>:	mulps  %xmm14,%xmm15
   0x0000000000001aa3 <+6819>:	mulps  %xmm14,%xmm13
   0x0000000000001ab7 <+6839>:	mulps  %xmm14,%xmm12

567	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_14]), C_row[i]);
   0x0000000000001a79 <+6777>:	addps  0xc0(%r9),%xmm0
   0x0000000000001a93 <+6803>:	addps  0xd0(%r9),%xmm15
   0x0000000000001aa7 <+6823>:	addps  0xe0(%r9),%xmm13
   0x0000000000001abb <+6843>:	addps  0xf0(%r9),%xmm12

568	        _mm_store_ps(&C[i*4+C_OFFSET_14], C_row[i]);
   0x0000000000001a81 <+6785>:	movaps %xmm0,0xc0(%r9)
   0x0000000000001a89 <+6793>:	movss  0x10(%rsp),%xmm0
   0x0000000000001a9b <+6811>:	movaps %xmm15,0xd0(%r9)
   0x0000000000001aaf <+6831>:	movaps %xmm13,0xe0(%r9)
   0x0000000000001ac3 <+6851>:	movaps %xmm12,0xf0(%r9)

569	      }
570	    }
571
572	    /* Reset C(2,1) matrix accumulators */
573	    C_row[0] = _mm_setzero_ps();
574	    C_row[1] = _mm_setzero_ps();
575	    C_row[2] = _mm_setzero_ps();
576	    C_row[3] = _mm_setzero_ps();
577
578	    if (norm_product[4][0] &&
   0x0000000000001acb <+6859>:	movss  0x148(%rsp),%xmm12
   0x0000000000001ad5 <+6869>:	comiss %xmm1,%xmm12
   0x0000000000001ad9 <+6873>:	jb     0x2032 <stream_kernel+8242>

579	        norm_product[5][4] &&
   0x0000000000001adf <+6879>:	movss  0x130(%rsp),%xmm12
   0x0000000000001ae9 <+6889>:	comiss %xmm1,%xmm12
   0x0000000000001aed <+6893>:	jb     0x2032 <stream_kernel+8242>

580	        norm_product[6][8] &&
   0x0000000000001af3 <+6899>:	movss  0x120(%rsp),%xmm12
   0x0000000000001afd <+6909>:	comiss %xmm1,%xmm12
   0x0000000000001b01 <+6913>:	jb     0x2032 <stream_kernel+8242>

581	        norm_product[7][12])
   0x0000000000001b07 <+6919>:	movss  0x118(%rsp),%xmm12
   0x0000000000001b11 <+6929>:	comiss %xmm1,%xmm12
   0x0000000000001b15 <+6933>:	jb     0x2032 <stream_kernel+8242>

582	    {
583	      /* A(2,1)*B(1,1) = C(2,1). */
584	      for (i = 0; i < 4; i++)
585	      {
586	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
587	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
588	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b1b <+6939>:	movaps 0x400(%r8),%xmm15
   0x0000000000001b2b <+6955>:	movaps 0x440(%r8),%xmm12
   0x0000000000001b3b <+6971>:	mulps  (%rcx),%xmm15
   0x0000000000001b44 <+6980>:	mulps  (%rcx),%xmm12
   0x0000000000001b74 <+7028>:	movaps 0x480(%r8),%xmm12
   0x0000000000001b7c <+7036>:	mulps  (%rcx),%xmm12
   0x0000000000001bea <+7146>:	movaps 0x4c0(%r8),%xmm14
   0x0000000000001bf2 <+7154>:	mulps  (%rcx),%xmm14

589
590	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
591	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
592	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b23 <+6947>:	movaps 0x410(%r8),%xmm13
   0x0000000000001b33 <+6963>:	movaps 0x450(%r8),%xmm14
   0x0000000000001b3f <+6975>:	mulps  0x10(%rcx),%xmm13
   0x0000000000001b48 <+6984>:	mulps  0x10(%rcx),%xmm14
   0x0000000000001b4d <+6989>:	movss  %xmm0,0x10(%rsp)
   0x0000000000001b5b <+7003>:	addps  %xmm15,%xmm13
   0x0000000000001b70 <+7024>:	addps  %xmm12,%xmm14
   0x0000000000001bb7 <+7095>:	movaps 0x490(%r8),%xmm13
   0x0000000000001bbf <+7103>:	mulps  0x10(%rcx),%xmm13
   0x0000000000001bc4 <+7108>:	addps  %xmm12,%xmm13
   0x0000000000001bc8 <+7112>:	movaps 0x4d0(%r8),%xmm12
   0x0000000000001bd0 <+7120>:	mulps  0x10(%rcx),%xmm12
   0x0000000000001bf6 <+7158>:	addps  %xmm14,%xmm12

593
594	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
595	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
596	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b5f <+7007>:	movaps 0x420(%r8),%xmm15
   0x0000000000001b67 <+7015>:	mulps  0x20(%rcx),%xmm15
   0x0000000000001b80 <+7040>:	addps  %xmm13,%xmm15
   0x0000000000001b84 <+7044>:	movaps 0x460(%r8),%xmm13
   0x0000000000001b8c <+7052>:	mulps  0x20(%rcx),%xmm13
   0x0000000000001ba2 <+7074>:	addps  %xmm14,%xmm13
   0x0000000000001ba6 <+7078>:	movaps 0x4a0(%r8),%xmm14
   0x0000000000001bae <+7086>:	mulps  0x20(%rcx),%xmm14
   0x0000000000001bd5 <+7125>:	addps  %xmm13,%xmm14
   0x0000000000001bfa <+7162>:	movaps 0x4e0(%r8),%xmm14
   0x0000000000001c02 <+7170>:	mulps  0x20(%rcx),%xmm14
   0x0000000000001c07 <+7175>:	addps  %xmm12,%xmm14

597
598	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
599	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
600	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b53 <+6995>:	movaps 0x430(%r8),%xmm0
   0x0000000000001b6c <+7020>:	mulps  0x30(%rcx),%xmm0
   0x0000000000001b91 <+7057>:	addps  %xmm15,%xmm0
   0x0000000000001b95 <+7061>:	movaps 0x470(%r8),%xmm15
   0x0000000000001b9d <+7069>:	mulps  0x30(%rcx),%xmm15
   0x0000000000001bb3 <+7091>:	addps  %xmm13,%xmm15
   0x0000000000001bd9 <+7129>:	movaps 0x4b0(%r8),%xmm13
   0x0000000000001be1 <+7137>:	mulps  0x30(%rcx),%xmm13
   0x0000000000001be6 <+7142>:	addps  %xmm14,%xmm13
   0x0000000000001c0b <+7179>:	movaps 0x4f0(%r8),%xmm12
   0x0000000000001c13 <+7187>:	mulps  0x30(%rcx),%xmm12
   0x0000000000001c18 <+7192>:	addps  %xmm14,%xmm12

601	      }
602
603	      /* A(2,2)*B(2,1) = C(2,1). */
604	      for (i = 0; i < 4; i++)
605	      {
606	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
607	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
608	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c1c <+7196>:	movaps 0x500(%r8),%xmm14
   0x0000000000001c24 <+7204>:	mulps  0x100(%rcx),%xmm14
   0x0000000000001c2c <+7212>:	addps  %xmm0,%xmm14
   0x0000000000001c6a <+7274>:	movaps 0x540(%r8),%xmm14
   0x0000000000001c72 <+7282>:	mulps  0x100(%rcx),%xmm14
   0x0000000000001c7a <+7290>:	addps  %xmm15,%xmm14
   0x0000000000001cba <+7354>:	movaps 0x580(%r8),%xmm14
   0x0000000000001cc2 <+7362>:	mulps  0x100(%rcx),%xmm14
   0x0000000000001cca <+7370>:	addps  %xmm13,%xmm14
   0x0000000000001d0a <+7434>:	movaps 0x5c0(%r8),%xmm14
   0x0000000000001d12 <+7442>:	mulps  0x100(%rcx),%xmm14
   0x0000000000001d1a <+7450>:	addps  %xmm12,%xmm14

609
610	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
611	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
612	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c30 <+7216>:	movaps 0x510(%r8),%xmm0
   0x0000000000001c38 <+7224>:	mulps  0x110(%rcx),%xmm0
   0x0000000000001c3f <+7231>:	addps  %xmm14,%xmm0
   0x0000000000001c7e <+7294>:	movaps 0x550(%r8),%xmm15
   0x0000000000001c86 <+7302>:	mulps  0x110(%rcx),%xmm15
   0x0000000000001c8e <+7310>:	addps  %xmm14,%xmm15
   0x0000000000001cce <+7374>:	movaps 0x590(%r8),%xmm13
   0x0000000000001cd6 <+7382>:	mulps  0x110(%rcx),%xmm13
   0x0000000000001cde <+7390>:	addps  %xmm14,%xmm13
   0x0000000000001d1e <+7454>:	movaps 0x5d0(%r8),%xmm12
   0x0000000000001d26 <+7462>:	mulps  0x110(%rcx),%xmm12
   0x0000000000001d2e <+7470>:	addps  %xmm14,%xmm12

613
614	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
615	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
616	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c43 <+7235>:	movaps 0x520(%r8),%xmm14
   0x0000000000001c4b <+7243>:	mulps  0x120(%rcx),%xmm14
   0x0000000000001c53 <+7251>:	addps  %xmm0,%xmm14
   0x0000000000001c92 <+7314>:	movaps 0x560(%r8),%xmm14
   0x0000000000001c9a <+7322>:	mulps  0x120(%rcx),%xmm14
   0x0000000000001ca2 <+7330>:	addps  %xmm15,%xmm14
   0x0000000000001ce2 <+7394>:	movaps 0x5a0(%r8),%xmm14
   0x0000000000001cea <+7402>:	mulps  0x120(%rcx),%xmm14
   0x0000000000001cf2 <+7410>:	addps  %xmm13,%xmm14
   0x0000000000001d32 <+7474>:	movaps 0x5e0(%r8),%xmm14
   0x0000000000001d3a <+7482>:	mulps  0x120(%rcx),%xmm14
   0x0000000000001d42 <+7490>:	addps  %xmm12,%xmm14

617
618	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
619	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
620	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c57 <+7255>:	movaps 0x530(%r8),%xmm0
   0x0000000000001c5f <+7263>:	mulps  0x130(%rcx),%xmm0
   0x0000000000001c66 <+7270>:	addps  %xmm14,%xmm0
   0x0000000000001ca6 <+7334>:	movaps 0x570(%r8),%xmm15
   0x0000000000001cae <+7342>:	mulps  0x130(%rcx),%xmm15
   0x0000000000001cb6 <+7350>:	addps  %xmm14,%xmm15
   0x0000000000001cf6 <+7414>:	movaps 0x5b0(%r8),%xmm13
   0x0000000000001cfe <+7422>:	mulps  0x130(%rcx),%xmm13
   0x0000000000001d06 <+7430>:	addps  %xmm14,%xmm13
   0x0000000000001d46 <+7494>:	movaps 0x5f0(%r8),%xmm12
   0x0000000000001d4e <+7502>:	mulps  0x130(%rcx),%xmm12
   0x0000000000001d56 <+7510>:	addps  %xmm14,%xmm12

621	      }
622
623	      /* A(2,3)*B(3,1) = C(2,1). */
624	      for (i = 0; i < 4; i++)
625	      {
626	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
627	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
628	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001d5a <+7514>:	movaps 0x600(%r8),%xmm14
   0x0000000000001d62 <+7522>:	mulps  0x200(%rcx),%xmm14
   0x0000000000001d6a <+7530>:	addps  %xmm0,%xmm14
   0x0000000000001da8 <+7592>:	movaps 0x640(%r8),%xmm14
   0x0000000000001db0 <+7600>:	mulps  0x200(%rcx),%xmm14
   0x0000000000001db8 <+7608>:	addps  %xmm15,%xmm14
   0x0000000000001df8 <+7672>:	movaps 0x680(%r8),%xmm14
   0x0000000000001e00 <+7680>:	mulps  0x200(%rcx),%xmm14
   0x0000000000001e08 <+7688>:	addps  %xmm13,%xmm14
   0x0000000000001e48 <+7752>:	movaps 0x6c0(%r8),%xmm14
   0x0000000000001e50 <+7760>:	mulps  0x200(%rcx),%xmm14
   0x0000000000001e58 <+7768>:	addps  %xmm12,%xmm14

629
630	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
631	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
632	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001d6e <+7534>:	movaps 0x610(%r8),%xmm0
   0x0000000000001d76 <+7542>:	mulps  0x210(%rcx),%xmm0
   0x0000000000001d7d <+7549>:	addps  %xmm14,%xmm0
   0x0000000000001dbc <+7612>:	movaps 0x650(%r8),%xmm15
   0x0000000000001dc4 <+7620>:	mulps  0x210(%rcx),%xmm15
   0x0000000000001dcc <+7628>:	addps  %xmm14,%xmm15
   0x0000000000001e0c <+7692>:	movaps 0x690(%r8),%xmm13
   0x0000000000001e14 <+7700>:	mulps  0x210(%rcx),%xmm13
   0x0000000000001e1c <+7708>:	addps  %xmm14,%xmm13
   0x0000000000001e5c <+7772>:	movaps 0x6d0(%r8),%xmm12
   0x0000000000001e64 <+7780>:	mulps  0x210(%rcx),%xmm12
   0x0000000000001e6c <+7788>:	addps  %xmm14,%xmm12

633
634	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
635	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
636	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001d81 <+7553>:	movaps 0x620(%r8),%xmm14
   0x0000000000001d89 <+7561>:	mulps  0x220(%rcx),%xmm14
   0x0000000000001d91 <+7569>:	addps  %xmm0,%xmm14
   0x0000000000001dd0 <+7632>:	movaps 0x660(%r8),%xmm14
   0x0000000000001dd8 <+7640>:	mulps  0x220(%rcx),%xmm14
   0x0000000000001de0 <+7648>:	addps  %xmm15,%xmm14
   0x0000000000001e20 <+7712>:	movaps 0x6a0(%r8),%xmm14
   0x0000000000001e28 <+7720>:	mulps  0x220(%rcx),%xmm14
   0x0000000000001e30 <+7728>:	addps  %xmm13,%xmm14
   0x0000000000001e70 <+7792>:	movaps 0x6e0(%r8),%xmm14
   0x0000000000001e78 <+7800>:	mulps  0x220(%rcx),%xmm14
   0x0000000000001e80 <+7808>:	addps  %xmm12,%xmm14

637
638	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
639	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
640	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001d95 <+7573>:	movaps 0x630(%r8),%xmm0
   0x0000000000001d9d <+7581>:	mulps  0x230(%rcx),%xmm0
   0x0000000000001da4 <+7588>:	addps  %xmm14,%xmm0
   0x0000000000001de4 <+7652>:	movaps 0x670(%r8),%xmm15
   0x0000000000001dec <+7660>:	mulps  0x230(%rcx),%xmm15
   0x0000000000001df4 <+7668>:	addps  %xmm14,%xmm15
   0x0000000000001e34 <+7732>:	movaps 0x6b0(%r8),%xmm13
   0x0000000000001e3c <+7740>:	mulps  0x230(%rcx),%xmm13
   0x0000000000001e44 <+7748>:	addps  %xmm14,%xmm13
   0x0000000000001e84 <+7812>:	movaps 0x6f0(%r8),%xmm12
   0x0000000000001e8c <+7820>:	mulps  0x230(%rcx),%xmm12
   0x0000000000001e94 <+7828>:	addps  %xmm14,%xmm12

641	      }
642
643	      /* A(2,4)*B(4,1) = C(2,1). */
644	      for (i = 0; i < 4; i++)
645	      {
646	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
647	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
648	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001e98 <+7832>:	movaps 0x700(%r8),%xmm14
   0x0000000000001ea0 <+7840>:	mulps  0x300(%rcx),%xmm14
   0x0000000000001ea8 <+7848>:	addps  %xmm0,%xmm14
   0x0000000000001ee6 <+7910>:	movaps 0x740(%r8),%xmm14
   0x0000000000001eee <+7918>:	mulps  0x300(%rcx),%xmm14
   0x0000000000001ef6 <+7926>:	addps  %xmm15,%xmm14
   0x0000000000001f36 <+7990>:	movaps 0x780(%r8),%xmm14
   0x0000000000001f3e <+7998>:	mulps  0x300(%rcx),%xmm14
   0x0000000000001f46 <+8006>:	addps  %xmm13,%xmm14
   0x0000000000001f86 <+8070>:	movaps 0x7c0(%r8),%xmm14
   0x0000000000001f8e <+8078>:	mulps  0x300(%rcx),%xmm14
   0x0000000000001f96 <+8086>:	addps  %xmm12,%xmm14

649
650	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
651	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
652	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001eac <+7852>:	movaps 0x710(%r8),%xmm0
   0x0000000000001eb4 <+7860>:	mulps  0x310(%rcx),%xmm0
   0x0000000000001ebb <+7867>:	addps  %xmm14,%xmm0
   0x0000000000001efa <+7930>:	movaps 0x750(%r8),%xmm15
   0x0000000000001f02 <+7938>:	mulps  0x310(%rcx),%xmm15
   0x0000000000001f0a <+7946>:	addps  %xmm14,%xmm15
   0x0000000000001f4a <+8010>:	movaps 0x790(%r8),%xmm13
   0x0000000000001f52 <+8018>:	mulps  0x310(%rcx),%xmm13
   0x0000000000001f5a <+8026>:	addps  %xmm14,%xmm13
   0x0000000000001f9a <+8090>:	movaps 0x7d0(%r8),%xmm12
   0x0000000000001fa2 <+8098>:	mulps  0x310(%rcx),%xmm12
   0x0000000000001faa <+8106>:	addps  %xmm14,%xmm12

653
654	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
655	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
656	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001ebf <+7871>:	movaps 0x720(%r8),%xmm14
   0x0000000000001ec7 <+7879>:	mulps  0x320(%rcx),%xmm14
   0x0000000000001ecf <+7887>:	addps  %xmm0,%xmm14
   0x0000000000001f0e <+7950>:	movaps 0x760(%r8),%xmm14
   0x0000000000001f16 <+7958>:	mulps  0x320(%rcx),%xmm14
   0x0000000000001f1e <+7966>:	addps  %xmm15,%xmm14
   0x0000000000001f5e <+8030>:	movaps 0x7a0(%r8),%xmm14
   0x0000000000001f66 <+8038>:	mulps  0x320(%rcx),%xmm14
   0x0000000000001f6e <+8046>:	addps  %xmm13,%xmm14
   0x0000000000001fae <+8110>:	movaps 0x7e0(%r8),%xmm14
   0x0000000000001fb6 <+8118>:	mulps  0x320(%rcx),%xmm14
   0x0000000000001fbe <+8126>:	addps  %xmm12,%xmm14

657
658	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
659	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
660	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001ed3 <+7891>:	movaps 0x730(%r8),%xmm0
   0x0000000000001edb <+7899>:	mulps  0x330(%rcx),%xmm0
   0x0000000000001ee2 <+7906>:	addps  %xmm14,%xmm0
   0x0000000000001f22 <+7970>:	movaps 0x770(%r8),%xmm15
   0x0000000000001f2a <+7978>:	mulps  0x330(%rcx),%xmm15
   0x0000000000001f32 <+7986>:	addps  %xmm14,%xmm15
   0x0000000000001f72 <+8050>:	movaps 0x7b0(%r8),%xmm13
   0x0000000000001f7a <+8058>:	mulps  0x330(%rcx),%xmm13
   0x0000000000001f82 <+8066>:	addps  %xmm14,%xmm13
   0x0000000000001fc2 <+8130>:	movaps 0x7f0(%r8),%xmm12
   0x0000000000001fca <+8138>:	mulps  0x330(%rcx),%xmm12
   0x0000000000001fd2 <+8146>:	addps  %xmm14,%xmm12

661	      }
662
663	      /* Store C(2,1) block. */
664	      for (i = 0; i < 4; i++)
665	      {
666	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000001fd6 <+8150>:	movaps 0x20(%rsp),%xmm14
   0x0000000000001fdc <+8156>:	mulps  %xmm14,%xmm0
   0x0000000000001ff6 <+8182>:	mulps  %xmm14,%xmm15
   0x000000000000200a <+8202>:	mulps  %xmm14,%xmm13
   0x000000000000201e <+8222>:	mulps  %xmm14,%xmm12

667	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_21]), C_row[i]);
   0x0000000000001fe0 <+8160>:	addps  0x100(%r9),%xmm0
   0x0000000000001ffa <+8186>:	addps  0x110(%r9),%xmm15
   0x000000000000200e <+8206>:	addps  0x120(%r9),%xmm13
   0x0000000000002022 <+8226>:	addps  0x130(%r9),%xmm12

668	        _mm_store_ps(&C[i*4+C_OFFSET_21], C_row[i]);
   0x0000000000001fe8 <+8168>:	movaps %xmm0,0x100(%r9)
   0x0000000000001ff0 <+8176>:	movss  0x10(%rsp),%xmm0
   0x0000000000002002 <+8194>:	movaps %xmm15,0x110(%r9)
   0x0000000000002016 <+8214>:	movaps %xmm13,0x120(%r9)
   0x000000000000202a <+8234>:	movaps %xmm12,0x130(%r9)

669	      }
670	    }
671
672	    /* Reset C(2,2) matrix accumulators */
673	    C_row[0] = _mm_setzero_ps();
674	    C_row[1] = _mm_setzero_ps();
675	    C_row[2] = _mm_setzero_ps();
676	    C_row[3] = _mm_setzero_ps();
677
678	    if (norm_product[4][1] &&
   0x0000000000002032 <+8242>:	comiss %xmm1,%xmm11
   0x0000000000002036 <+8246>:	jb     0x258e <stream_kernel+9614>

679	        norm_product[5][5] &&
   0x000000000000203c <+8252>:	movss  0x110(%rsp),%xmm11
   0x0000000000002046 <+8262>:	comiss %xmm1,%xmm11
   0x000000000000204a <+8266>:	jb     0x258e <stream_kernel+9614>

680	        norm_product[6][9] &&
   0x0000000000002050 <+8272>:	movss  0x108(%rsp),%xmm11
   0x000000000000205a <+8282>:	comiss %xmm1,%xmm11
   0x000000000000205e <+8286>:	jb     0x258e <stream_kernel+9614>

681	        norm_product[7][13])
   0x0000000000002064 <+8292>:	movss  0x100(%rsp),%xmm11
   0x000000000000206e <+8302>:	comiss %xmm1,%xmm11
   0x0000000000002072 <+8306>:	jb     0x258e <stream_kernel+9614>

682	    {
683	      /* A(2,1)*B(1,2) = C(2,2). */
684	      for (i = 0; i < 4; i++)
685	      {
686	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
687	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
688	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002078 <+8312>:	movaps 0x400(%r8),%xmm15
   0x0000000000002090 <+8336>:	movaps 0x440(%r8),%xmm12
   0x00000000000020a0 <+8352>:	mulps  0x40(%rcx),%xmm15
   0x00000000000020af <+8367>:	mulps  0x40(%rcx),%xmm12
   0x00000000000020ce <+8398>:	movaps 0x480(%r8),%xmm12
   0x00000000000020d6 <+8406>:	mulps  0x40(%rcx),%xmm12
   0x0000000000002145 <+8517>:	movaps 0x4c0(%r8),%xmm14
   0x000000000000214d <+8525>:	mulps  0x40(%rcx),%xmm14

689
690	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
691	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
692	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002080 <+8320>:	movaps 0x410(%r8),%xmm13
   0x0000000000002098 <+8344>:	movaps 0x450(%r8),%xmm14
   0x00000000000020a5 <+8357>:	mulps  0x50(%rcx),%xmm13
   0x00000000000020b4 <+8372>:	mulps  0x50(%rcx),%xmm14
   0x00000000000020b9 <+8377>:	addps  %xmm15,%xmm13
   0x00000000000020ca <+8394>:	addps  %xmm12,%xmm14
   0x0000000000002112 <+8466>:	movaps 0x490(%r8),%xmm13
   0x000000000000211a <+8474>:	mulps  0x50(%rcx),%xmm13
   0x000000000000211f <+8479>:	addps  %xmm12,%xmm13
   0x0000000000002123 <+8483>:	movaps 0x4d0(%r8),%xmm12
   0x000000000000212b <+8491>:	mulps  0x50(%rcx),%xmm12
   0x0000000000002152 <+8530>:	addps  %xmm14,%xmm12

693
694	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
695	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
696	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000020bd <+8381>:	movaps 0x420(%r8),%xmm15
   0x00000000000020c5 <+8389>:	mulps  0x60(%rcx),%xmm15
   0x00000000000020db <+8411>:	addps  %xmm13,%xmm15
   0x00000000000020df <+8415>:	movaps 0x460(%r8),%xmm13
   0x00000000000020e7 <+8423>:	mulps  0x60(%rcx),%xmm13
   0x00000000000020fd <+8445>:	addps  %xmm14,%xmm13
   0x0000000000002101 <+8449>:	movaps 0x4a0(%r8),%xmm14
   0x0000000000002109 <+8457>:	mulps  0x60(%rcx),%xmm14
   0x0000000000002130 <+8496>:	addps  %xmm13,%xmm14
   0x0000000000002156 <+8534>:	movaps 0x4e0(%r8),%xmm14
   0x000000000000215e <+8542>:	mulps  0x60(%rcx),%xmm14
   0x0000000000002163 <+8547>:	addps  %xmm12,%xmm14

697
698	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
699	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
700	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002088 <+8328>:	movaps 0x430(%r8),%xmm11
   0x00000000000020aa <+8362>:	mulps  0x70(%rcx),%xmm11
   0x00000000000020ec <+8428>:	addps  %xmm15,%xmm11
   0x00000000000020f0 <+8432>:	movaps 0x470(%r8),%xmm15
   0x00000000000020f8 <+8440>:	mulps  0x70(%rcx),%xmm15
   0x000000000000210e <+8462>:	addps  %xmm13,%xmm15
   0x0000000000002134 <+8500>:	movaps 0x4b0(%r8),%xmm13
   0x000000000000213c <+8508>:	mulps  0x70(%rcx),%xmm13
   0x0000000000002141 <+8513>:	addps  %xmm14,%xmm13
   0x0000000000002167 <+8551>:	movaps 0x4f0(%r8),%xmm12
   0x000000000000216f <+8559>:	mulps  0x70(%rcx),%xmm12
   0x0000000000002174 <+8564>:	addps  %xmm14,%xmm12

701	      }
702
703	      /* A(2,2)*B(2,2) = C(2,2). */
704	      for (i = 0; i < 4; i++)
705	      {
706	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
707	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
708	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002178 <+8568>:	movaps 0x500(%r8),%xmm14
   0x0000000000002180 <+8576>:	mulps  0x140(%rcx),%xmm14
   0x0000000000002188 <+8584>:	addps  %xmm11,%xmm14
   0x00000000000021c8 <+8648>:	movaps 0x540(%r8),%xmm14
   0x00000000000021d0 <+8656>:	mulps  0x140(%rcx),%xmm14
   0x00000000000021d8 <+8664>:	addps  %xmm15,%xmm14
   0x0000000000002218 <+8728>:	movaps 0x580(%r8),%xmm14
   0x0000000000002220 <+8736>:	mulps  0x140(%rcx),%xmm14
   0x0000000000002228 <+8744>:	addps  %xmm13,%xmm14
   0x0000000000002268 <+8808>:	movaps 0x5c0(%r8),%xmm14
   0x0000000000002270 <+8816>:	mulps  0x140(%rcx),%xmm14
   0x0000000000002278 <+8824>:	addps  %xmm12,%xmm14

709
710	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
711	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
712	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000218c <+8588>:	movaps 0x510(%r8),%xmm11
   0x0000000000002194 <+8596>:	mulps  0x150(%rcx),%xmm11
   0x000000000000219c <+8604>:	addps  %xmm14,%xmm11
   0x00000000000021dc <+8668>:	movaps 0x550(%r8),%xmm15
   0x00000000000021e4 <+8676>:	mulps  0x150(%rcx),%xmm15
   0x00000000000021ec <+8684>:	addps  %xmm14,%xmm15
   0x000000000000222c <+8748>:	movaps 0x590(%r8),%xmm13
   0x0000000000002234 <+8756>:	mulps  0x150(%rcx),%xmm13
   0x000000000000223c <+8764>:	addps  %xmm14,%xmm13
   0x000000000000227c <+8828>:	movaps 0x5d0(%r8),%xmm12
   0x0000000000002284 <+8836>:	mulps  0x150(%rcx),%xmm12
   0x000000000000228c <+8844>:	addps  %xmm14,%xmm12

713
714	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
715	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
716	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000021a0 <+8608>:	movaps 0x520(%r8),%xmm14
   0x00000000000021a8 <+8616>:	mulps  0x160(%rcx),%xmm14
   0x00000000000021b0 <+8624>:	addps  %xmm11,%xmm14
   0x00000000000021f0 <+8688>:	movaps 0x560(%r8),%xmm14
   0x00000000000021f8 <+8696>:	mulps  0x160(%rcx),%xmm14
   0x0000000000002200 <+8704>:	addps  %xmm15,%xmm14
   0x0000000000002240 <+8768>:	movaps 0x5a0(%r8),%xmm14
   0x0000000000002248 <+8776>:	mulps  0x160(%rcx),%xmm14
   0x0000000000002250 <+8784>:	addps  %xmm13,%xmm14
   0x0000000000002290 <+8848>:	movaps 0x5e0(%r8),%xmm14
   0x0000000000002298 <+8856>:	mulps  0x160(%rcx),%xmm14
   0x00000000000022a0 <+8864>:	addps  %xmm12,%xmm14

717
718	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
719	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
720	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000021b4 <+8628>:	movaps 0x530(%r8),%xmm11
   0x00000000000021bc <+8636>:	mulps  0x170(%rcx),%xmm11
   0x00000000000021c4 <+8644>:	addps  %xmm14,%xmm11
   0x0000000000002204 <+8708>:	movaps 0x570(%r8),%xmm15
   0x000000000000220c <+8716>:	mulps  0x170(%rcx),%xmm15
   0x0000000000002214 <+8724>:	addps  %xmm14,%xmm15
   0x0000000000002254 <+8788>:	movaps 0x5b0(%r8),%xmm13
   0x000000000000225c <+8796>:	mulps  0x170(%rcx),%xmm13
   0x0000000000002264 <+8804>:	addps  %xmm14,%xmm13
   0x00000000000022a4 <+8868>:	movaps 0x5f0(%r8),%xmm12
   0x00000000000022ac <+8876>:	mulps  0x170(%rcx),%xmm12
   0x00000000000022b4 <+8884>:	addps  %xmm14,%xmm12

721	      }
722
723	      /* A(2,3)*B(3,2) = C(2,2). */
724	      for (i = 0; i < 4; i++)
725	      {
726	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
727	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
728	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000022b8 <+8888>:	movaps 0x600(%r8),%xmm14
   0x00000000000022c0 <+8896>:	mulps  0x240(%rcx),%xmm14
   0x00000000000022c8 <+8904>:	addps  %xmm11,%xmm14
   0x0000000000002308 <+8968>:	movaps 0x640(%r8),%xmm14
   0x0000000000002310 <+8976>:	mulps  0x240(%rcx),%xmm14
   0x0000000000002318 <+8984>:	addps  %xmm15,%xmm14
   0x0000000000002358 <+9048>:	movaps 0x680(%r8),%xmm14
   0x0000000000002360 <+9056>:	mulps  0x240(%rcx),%xmm14
   0x0000000000002368 <+9064>:	addps  %xmm13,%xmm14
   0x00000000000023a8 <+9128>:	movaps 0x6c0(%r8),%xmm14
   0x00000000000023b0 <+9136>:	mulps  0x240(%rcx),%xmm14
   0x00000000000023b8 <+9144>:	addps  %xmm12,%xmm14

729
730	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
731	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
732	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000022cc <+8908>:	movaps 0x610(%r8),%xmm11
   0x00000000000022d4 <+8916>:	mulps  0x250(%rcx),%xmm11
   0x00000000000022dc <+8924>:	addps  %xmm14,%xmm11
   0x000000000000231c <+8988>:	movaps 0x650(%r8),%xmm15
   0x0000000000002324 <+8996>:	mulps  0x250(%rcx),%xmm15
   0x000000000000232c <+9004>:	addps  %xmm14,%xmm15
   0x000000000000236c <+9068>:	movaps 0x690(%r8),%xmm13
   0x0000000000002374 <+9076>:	mulps  0x250(%rcx),%xmm13
   0x000000000000237c <+9084>:	addps  %xmm14,%xmm13
   0x00000000000023bc <+9148>:	movaps 0x6d0(%r8),%xmm12
   0x00000000000023c4 <+9156>:	mulps  0x250(%rcx),%xmm12
   0x00000000000023cc <+9164>:	addps  %xmm14,%xmm12

733
734	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
735	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
736	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000022e0 <+8928>:	movaps 0x620(%r8),%xmm14
   0x00000000000022e8 <+8936>:	mulps  0x260(%rcx),%xmm14
   0x00000000000022f0 <+8944>:	addps  %xmm11,%xmm14
   0x0000000000002330 <+9008>:	movaps 0x660(%r8),%xmm14
   0x0000000000002338 <+9016>:	mulps  0x260(%rcx),%xmm14
   0x0000000000002340 <+9024>:	addps  %xmm15,%xmm14
   0x0000000000002380 <+9088>:	movaps 0x6a0(%r8),%xmm14
   0x0000000000002388 <+9096>:	mulps  0x260(%rcx),%xmm14
   0x0000000000002390 <+9104>:	addps  %xmm13,%xmm14
   0x00000000000023d0 <+9168>:	movaps 0x6e0(%r8),%xmm14
   0x00000000000023d8 <+9176>:	mulps  0x260(%rcx),%xmm14
   0x00000000000023e0 <+9184>:	addps  %xmm12,%xmm14

737
738	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
739	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
740	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000022f4 <+8948>:	movaps 0x630(%r8),%xmm11
   0x00000000000022fc <+8956>:	mulps  0x270(%rcx),%xmm11
   0x0000000000002304 <+8964>:	addps  %xmm14,%xmm11
   0x0000000000002344 <+9028>:	movaps 0x670(%r8),%xmm15
   0x000000000000234c <+9036>:	mulps  0x270(%rcx),%xmm15
   0x0000000000002354 <+9044>:	addps  %xmm14,%xmm15
   0x0000000000002394 <+9108>:	movaps 0x6b0(%r8),%xmm13
   0x000000000000239c <+9116>:	mulps  0x270(%rcx),%xmm13
   0x00000000000023a4 <+9124>:	addps  %xmm14,%xmm13
   0x00000000000023e4 <+9188>:	movaps 0x6f0(%r8),%xmm12
   0x00000000000023ec <+9196>:	mulps  0x270(%rcx),%xmm12
   0x00000000000023f4 <+9204>:	addps  %xmm14,%xmm12

741	      }
742
743	      /* A(2,4)*B(4,2) = C(2,2). */
744	      for (i = 0; i < 4; i++)
745	      {
746	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
747	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
748	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000023f8 <+9208>:	movaps 0x700(%r8),%xmm14
   0x0000000000002400 <+9216>:	mulps  0x340(%rcx),%xmm14
   0x0000000000002408 <+9224>:	addps  %xmm11,%xmm14
   0x0000000000002448 <+9288>:	movaps 0x740(%r8),%xmm14
   0x0000000000002450 <+9296>:	mulps  0x340(%rcx),%xmm14
   0x0000000000002458 <+9304>:	addps  %xmm15,%xmm14
   0x0000000000002498 <+9368>:	movaps 0x780(%r8),%xmm14
   0x00000000000024a0 <+9376>:	mulps  0x340(%rcx),%xmm14
   0x00000000000024a8 <+9384>:	addps  %xmm13,%xmm14
   0x00000000000024e8 <+9448>:	movaps 0x7c0(%r8),%xmm14
   0x00000000000024f0 <+9456>:	mulps  0x340(%rcx),%xmm14
   0x00000000000024f8 <+9464>:	addps  %xmm12,%xmm14

749
750	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
751	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
752	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000240c <+9228>:	movaps 0x710(%r8),%xmm11
   0x0000000000002414 <+9236>:	mulps  0x350(%rcx),%xmm11
   0x000000000000241c <+9244>:	addps  %xmm14,%xmm11
   0x000000000000245c <+9308>:	movaps 0x750(%r8),%xmm15
   0x0000000000002464 <+9316>:	mulps  0x350(%rcx),%xmm15
   0x000000000000246c <+9324>:	addps  %xmm14,%xmm15
   0x00000000000024ac <+9388>:	movaps 0x790(%r8),%xmm13
   0x00000000000024b4 <+9396>:	mulps  0x350(%rcx),%xmm13
   0x00000000000024bc <+9404>:	addps  %xmm14,%xmm13
   0x00000000000024fc <+9468>:	movaps 0x7d0(%r8),%xmm12
   0x0000000000002504 <+9476>:	mulps  0x350(%rcx),%xmm12
   0x000000000000250c <+9484>:	addps  %xmm14,%xmm12

753
754	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
755	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
756	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002420 <+9248>:	movaps 0x720(%r8),%xmm14
   0x0000000000002428 <+9256>:	mulps  0x360(%rcx),%xmm14
   0x0000000000002430 <+9264>:	addps  %xmm11,%xmm14
   0x0000000000002470 <+9328>:	movaps 0x760(%r8),%xmm14
   0x0000000000002478 <+9336>:	mulps  0x360(%rcx),%xmm14
   0x0000000000002480 <+9344>:	addps  %xmm15,%xmm14
   0x00000000000024c0 <+9408>:	movaps 0x7a0(%r8),%xmm14
   0x00000000000024c8 <+9416>:	mulps  0x360(%rcx),%xmm14
   0x00000000000024d0 <+9424>:	addps  %xmm13,%xmm14
   0x0000000000002510 <+9488>:	movaps 0x7e0(%r8),%xmm14
   0x0000000000002518 <+9496>:	mulps  0x360(%rcx),%xmm14
   0x0000000000002520 <+9504>:	addps  %xmm12,%xmm14

757
758	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
759	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
760	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002434 <+9268>:	movaps 0x730(%r8),%xmm11
   0x000000000000243c <+9276>:	mulps  0x370(%rcx),%xmm11
   0x0000000000002444 <+9284>:	addps  %xmm14,%xmm11
   0x0000000000002484 <+9348>:	movaps 0x770(%r8),%xmm15
   0x000000000000248c <+9356>:	mulps  0x370(%rcx),%xmm15
   0x0000000000002494 <+9364>:	addps  %xmm14,%xmm15
   0x00000000000024d4 <+9428>:	movaps 0x7b0(%r8),%xmm13
   0x00000000000024dc <+9436>:	mulps  0x370(%rcx),%xmm13
   0x00000000000024e4 <+9444>:	addps  %xmm14,%xmm13
   0x0000000000002524 <+9508>:	movaps 0x7f0(%r8),%xmm12
   0x000000000000252c <+9516>:	mulps  0x370(%rcx),%xmm12
   0x0000000000002534 <+9524>:	addps  %xmm14,%xmm12

761	      }
762
763	      /* Store C(2,2) block. */
764	      for (i = 0; i < 4; i++)
765	      {
766	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000002538 <+9528>:	movaps 0x20(%rsp),%xmm14
   0x000000000000253e <+9534>:	mulps  %xmm14,%xmm11
   0x0000000000002552 <+9554>:	mulps  %xmm14,%xmm15
   0x0000000000002566 <+9574>:	mulps  %xmm14,%xmm13
   0x000000000000257a <+9594>:	mulps  %xmm14,%xmm12

767	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_22]), C_row[i]);
   0x0000000000002542 <+9538>:	addps  0x140(%r9),%xmm11
   0x0000000000002556 <+9558>:	addps  0x150(%r9),%xmm15
   0x000000000000256a <+9578>:	addps  0x160(%r9),%xmm13
   0x000000000000257e <+9598>:	addps  0x170(%r9),%xmm12

768	        _mm_store_ps(&C[i*4+C_OFFSET_22], C_row[i]);
   0x000000000000254a <+9546>:	movaps %xmm11,0x140(%r9)
   0x000000000000255e <+9566>:	movaps %xmm15,0x150(%r9)
   0x0000000000002572 <+9586>:	movaps %xmm13,0x160(%r9)
   0x0000000000002586 <+9606>:	movaps %xmm12,0x170(%r9)

769	      }
770	    }
771
772	    /* Reset C(2,3) matrix accumulators */
773	    C_row[0] = _mm_setzero_ps();
774	    C_row[1] = _mm_setzero_ps();
775	    C_row[2] = _mm_setzero_ps();
776	    C_row[3] = _mm_setzero_ps();
777
778	    if (norm_product[4][2] &&
   0x000000000000258e <+9614>:	movss  0x178(%rsp),%xmm11
   0x0000000000002598 <+9624>:	comiss %xmm1,%xmm11
   0x000000000000259c <+9628>:	jb     0x2b24 <stream_kernel+11044>

779	        norm_product[5][6] &&
   0x00000000000025a2 <+9634>:	movss  0xf8(%rsp),%xmm11
   0x00000000000025ac <+9644>:	comiss %xmm1,%xmm11
   0x00000000000025b0 <+9648>:	jb     0x2b24 <stream_kernel+11044>

780	        norm_product[6][10] &&
   0x00000000000025b6 <+9654>:	movss  0xe8(%rsp),%xmm11
   0x00000000000025c0 <+9664>:	comiss %xmm1,%xmm11
   0x00000000000025c4 <+9668>:	jb     0x2b24 <stream_kernel+11044>

781	        norm_product[7][14])
   0x00000000000025ca <+9674>:	movss  0xe0(%rsp),%xmm11
   0x00000000000025d4 <+9684>:	comiss %xmm1,%xmm11
   0x00000000000025d8 <+9688>:	jb     0x2b24 <stream_kernel+11044>

782	    {
783	      /* A(2,1)*B(1,3) = C(2,3). */
784	      for (i = 0; i < 4; i++)
785	      {
786	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
787	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
788	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000025de <+9694>:	movaps 0x400(%r8),%xmm15
   0x00000000000025f6 <+9718>:	movaps 0x440(%r8),%xmm12
   0x0000000000002606 <+9734>:	mulps  0x80(%rcx),%xmm15
   0x000000000000261e <+9758>:	mulps  0x80(%rcx),%xmm12
   0x0000000000002646 <+9798>:	movaps 0x480(%r8),%xmm12
   0x000000000000264e <+9806>:	mulps  0x80(%rcx),%xmm12
   0x00000000000026d2 <+9938>:	movaps 0x4c0(%r8),%xmm14
   0x00000000000026da <+9946>:	mulps  0x80(%rcx),%xmm14

789
790	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
791	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
792	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000025e6 <+9702>:	movaps 0x410(%r8),%xmm13
   0x00000000000025fe <+9726>:	movaps 0x450(%r8),%xmm14
   0x000000000000260e <+9742>:	mulps  0x90(%rcx),%xmm13
   0x0000000000002626 <+9766>:	mulps  0x90(%rcx),%xmm14
   0x000000000000262e <+9774>:	addps  %xmm15,%xmm13
   0x0000000000002642 <+9794>:	addps  %xmm12,%xmm14
   0x0000000000002696 <+9878>:	movaps 0x490(%r8),%xmm13
   0x000000000000269e <+9886>:	mulps  0x90(%rcx),%xmm13
   0x00000000000026a6 <+9894>:	addps  %xmm12,%xmm13
   0x00000000000026aa <+9898>:	movaps 0x4d0(%r8),%xmm12
   0x00000000000026b2 <+9906>:	mulps  0x90(%rcx),%xmm12
   0x00000000000026e2 <+9954>:	addps  %xmm14,%xmm12

793
794	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
795	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
796	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002632 <+9778>:	movaps 0x420(%r8),%xmm15
   0x000000000000263a <+9786>:	mulps  0xa0(%rcx),%xmm15
   0x0000000000002656 <+9814>:	addps  %xmm13,%xmm15
   0x000000000000265a <+9818>:	movaps 0x460(%r8),%xmm13
   0x0000000000002662 <+9826>:	mulps  0xa0(%rcx),%xmm13
   0x000000000000267e <+9854>:	addps  %xmm14,%xmm13
   0x0000000000002682 <+9858>:	movaps 0x4a0(%r8),%xmm14
   0x000000000000268a <+9866>:	mulps  0xa0(%rcx),%xmm14
   0x00000000000026ba <+9914>:	addps  %xmm13,%xmm14
   0x00000000000026e6 <+9958>:	movaps 0x4e0(%r8),%xmm14
   0x00000000000026ee <+9966>:	mulps  0xa0(%rcx),%xmm14
   0x00000000000026f6 <+9974>:	addps  %xmm12,%xmm14

797
798	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
799	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
800	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000025ee <+9710>:	movaps 0x430(%r8),%xmm11
   0x0000000000002616 <+9750>:	mulps  0xb0(%rcx),%xmm11
   0x000000000000266a <+9834>:	addps  %xmm15,%xmm11
   0x000000000000266e <+9838>:	movaps 0x470(%r8),%xmm15
   0x0000000000002676 <+9846>:	mulps  0xb0(%rcx),%xmm15
   0x0000000000002692 <+9874>:	addps  %xmm13,%xmm15
   0x00000000000026be <+9918>:	movaps 0x4b0(%r8),%xmm13
   0x00000000000026c6 <+9926>:	mulps  0xb0(%rcx),%xmm13
   0x00000000000026ce <+9934>:	addps  %xmm14,%xmm13
   0x00000000000026fa <+9978>:	movaps 0x4f0(%r8),%xmm12
   0x0000000000002702 <+9986>:	mulps  0xb0(%rcx),%xmm12
   0x000000000000270a <+9994>:	addps  %xmm14,%xmm12

801	      }
802
803	      /* A(2,2)*B(2,3) = C(2,3). */
804	      for (i = 0; i < 4; i++)
805	      {
806	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
807	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
808	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000270e <+9998>:	movaps 0x500(%r8),%xmm14
   0x0000000000002716 <+10006>:	mulps  0x180(%rcx),%xmm14
   0x000000000000271e <+10014>:	addps  %xmm11,%xmm14
   0x000000000000275e <+10078>:	movaps 0x540(%r8),%xmm14
   0x0000000000002766 <+10086>:	mulps  0x180(%rcx),%xmm14
   0x000000000000276e <+10094>:	addps  %xmm15,%xmm14
   0x00000000000027ae <+10158>:	movaps 0x580(%r8),%xmm14
   0x00000000000027b6 <+10166>:	mulps  0x180(%rcx),%xmm14
   0x00000000000027be <+10174>:	addps  %xmm13,%xmm14
   0x00000000000027fe <+10238>:	movaps 0x5c0(%r8),%xmm14
   0x0000000000002806 <+10246>:	mulps  0x180(%rcx),%xmm14
   0x000000000000280e <+10254>:	addps  %xmm12,%xmm14

809
810	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
811	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
812	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002722 <+10018>:	movaps 0x510(%r8),%xmm11
   0x000000000000272a <+10026>:	mulps  0x190(%rcx),%xmm11
   0x0000000000002732 <+10034>:	addps  %xmm14,%xmm11
   0x0000000000002772 <+10098>:	movaps 0x550(%r8),%xmm15
   0x000000000000277a <+10106>:	mulps  0x190(%rcx),%xmm15
   0x0000000000002782 <+10114>:	addps  %xmm14,%xmm15
   0x00000000000027c2 <+10178>:	movaps 0x590(%r8),%xmm13
   0x00000000000027ca <+10186>:	mulps  0x190(%rcx),%xmm13
   0x00000000000027d2 <+10194>:	addps  %xmm14,%xmm13
   0x0000000000002812 <+10258>:	movaps 0x5d0(%r8),%xmm12
   0x000000000000281a <+10266>:	mulps  0x190(%rcx),%xmm12
   0x0000000000002822 <+10274>:	addps  %xmm14,%xmm12

813
814	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
815	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
816	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002736 <+10038>:	movaps 0x520(%r8),%xmm14
   0x000000000000273e <+10046>:	mulps  0x1a0(%rcx),%xmm14
   0x0000000000002746 <+10054>:	addps  %xmm11,%xmm14
   0x0000000000002786 <+10118>:	movaps 0x560(%r8),%xmm14
   0x000000000000278e <+10126>:	mulps  0x1a0(%rcx),%xmm14
   0x0000000000002796 <+10134>:	addps  %xmm15,%xmm14
   0x00000000000027d6 <+10198>:	movaps 0x5a0(%r8),%xmm14
   0x00000000000027de <+10206>:	mulps  0x1a0(%rcx),%xmm14
   0x00000000000027e6 <+10214>:	addps  %xmm13,%xmm14
   0x0000000000002826 <+10278>:	movaps 0x5e0(%r8),%xmm14
   0x000000000000282e <+10286>:	mulps  0x1a0(%rcx),%xmm14
   0x0000000000002836 <+10294>:	addps  %xmm12,%xmm14

817
818	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
819	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
820	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000274a <+10058>:	movaps 0x530(%r8),%xmm11
   0x0000000000002752 <+10066>:	mulps  0x1b0(%rcx),%xmm11
   0x000000000000275a <+10074>:	addps  %xmm14,%xmm11
   0x000000000000279a <+10138>:	movaps 0x570(%r8),%xmm15
   0x00000000000027a2 <+10146>:	mulps  0x1b0(%rcx),%xmm15
   0x00000000000027aa <+10154>:	addps  %xmm14,%xmm15
   0x00000000000027ea <+10218>:	movaps 0x5b0(%r8),%xmm13
   0x00000000000027f2 <+10226>:	mulps  0x1b0(%rcx),%xmm13
   0x00000000000027fa <+10234>:	addps  %xmm14,%xmm13
   0x000000000000283a <+10298>:	movaps 0x5f0(%r8),%xmm12
   0x0000000000002842 <+10306>:	mulps  0x1b0(%rcx),%xmm12
   0x000000000000284a <+10314>:	addps  %xmm14,%xmm12

821	      }
822
823	      /* A(2,3)*B(3,3) = C(2,3). */
824	      for (i = 0; i < 4; i++)
825	      {
826	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
827	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
828	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000284e <+10318>:	movaps 0x600(%r8),%xmm14
   0x0000000000002856 <+10326>:	mulps  0x280(%rcx),%xmm14
   0x000000000000285e <+10334>:	addps  %xmm11,%xmm14
   0x000000000000289e <+10398>:	movaps 0x640(%r8),%xmm14
   0x00000000000028a6 <+10406>:	mulps  0x280(%rcx),%xmm14
   0x00000000000028ae <+10414>:	addps  %xmm15,%xmm14
   0x00000000000028ee <+10478>:	movaps 0x680(%r8),%xmm14
   0x00000000000028f6 <+10486>:	mulps  0x280(%rcx),%xmm14
   0x00000000000028fe <+10494>:	addps  %xmm13,%xmm14
   0x000000000000293e <+10558>:	movaps 0x6c0(%r8),%xmm14
   0x0000000000002946 <+10566>:	mulps  0x280(%rcx),%xmm14
   0x000000000000294e <+10574>:	addps  %xmm12,%xmm14

829
830	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
831	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
832	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002862 <+10338>:	movaps 0x610(%r8),%xmm11
   0x000000000000286a <+10346>:	mulps  0x290(%rcx),%xmm11
   0x0000000000002872 <+10354>:	addps  %xmm14,%xmm11
   0x00000000000028b2 <+10418>:	movaps 0x650(%r8),%xmm15
   0x00000000000028ba <+10426>:	mulps  0x290(%rcx),%xmm15
   0x00000000000028c2 <+10434>:	addps  %xmm14,%xmm15
   0x0000000000002902 <+10498>:	movaps 0x690(%r8),%xmm13
   0x000000000000290a <+10506>:	mulps  0x290(%rcx),%xmm13
   0x0000000000002912 <+10514>:	addps  %xmm14,%xmm13
   0x0000000000002952 <+10578>:	movaps 0x6d0(%r8),%xmm12
   0x000000000000295a <+10586>:	mulps  0x290(%rcx),%xmm12
   0x0000000000002962 <+10594>:	addps  %xmm14,%xmm12

833
834	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
835	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
836	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002876 <+10358>:	movaps 0x620(%r8),%xmm14
   0x000000000000287e <+10366>:	mulps  0x2a0(%rcx),%xmm14
   0x0000000000002886 <+10374>:	addps  %xmm11,%xmm14
   0x00000000000028c6 <+10438>:	movaps 0x660(%r8),%xmm14
   0x00000000000028ce <+10446>:	mulps  0x2a0(%rcx),%xmm14
   0x00000000000028d6 <+10454>:	addps  %xmm15,%xmm14
   0x0000000000002916 <+10518>:	movaps 0x6a0(%r8),%xmm14
   0x000000000000291e <+10526>:	mulps  0x2a0(%rcx),%xmm14
   0x0000000000002926 <+10534>:	addps  %xmm13,%xmm14
   0x0000000000002966 <+10598>:	movaps 0x6e0(%r8),%xmm14
   0x000000000000296e <+10606>:	mulps  0x2a0(%rcx),%xmm14
   0x0000000000002976 <+10614>:	addps  %xmm12,%xmm14

837
838	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
839	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
840	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000288a <+10378>:	movaps 0x630(%r8),%xmm11
   0x0000000000002892 <+10386>:	mulps  0x2b0(%rcx),%xmm11
   0x000000000000289a <+10394>:	addps  %xmm14,%xmm11
   0x00000000000028da <+10458>:	movaps 0x670(%r8),%xmm15
   0x00000000000028e2 <+10466>:	mulps  0x2b0(%rcx),%xmm15
   0x00000000000028ea <+10474>:	addps  %xmm14,%xmm15
   0x000000000000292a <+10538>:	movaps 0x6b0(%r8),%xmm13
   0x0000000000002932 <+10546>:	mulps  0x2b0(%rcx),%xmm13
   0x000000000000293a <+10554>:	addps  %xmm14,%xmm13
   0x000000000000297a <+10618>:	movaps 0x6f0(%r8),%xmm12
   0x0000000000002982 <+10626>:	mulps  0x2b0(%rcx),%xmm12
   0x000000000000298a <+10634>:	addps  %xmm14,%xmm12

841	      }
842
843	      /* A(2,4)*B(4,3) = C(2,3). */
844	      for (i = 0; i < 4; i++)
845	      {
846	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
847	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
848	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000298e <+10638>:	movaps 0x700(%r8),%xmm14
   0x0000000000002996 <+10646>:	mulps  0x380(%rcx),%xmm14
   0x000000000000299e <+10654>:	addps  %xmm11,%xmm14
   0x00000000000029de <+10718>:	movaps 0x740(%r8),%xmm14
   0x00000000000029e6 <+10726>:	mulps  0x380(%rcx),%xmm14
   0x00000000000029ee <+10734>:	addps  %xmm15,%xmm14
   0x0000000000002a2e <+10798>:	movaps 0x780(%r8),%xmm14
   0x0000000000002a36 <+10806>:	mulps  0x380(%rcx),%xmm14
   0x0000000000002a3e <+10814>:	addps  %xmm13,%xmm14
   0x0000000000002a7e <+10878>:	movaps 0x7c0(%r8),%xmm14
   0x0000000000002a86 <+10886>:	mulps  0x380(%rcx),%xmm14
   0x0000000000002a8e <+10894>:	addps  %xmm12,%xmm14

849
850	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
851	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
852	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000029a2 <+10658>:	movaps 0x710(%r8),%xmm11
   0x00000000000029aa <+10666>:	mulps  0x390(%rcx),%xmm11
   0x00000000000029b2 <+10674>:	addps  %xmm14,%xmm11
   0x00000000000029f2 <+10738>:	movaps 0x750(%r8),%xmm15
   0x00000000000029fa <+10746>:	mulps  0x390(%rcx),%xmm15
   0x0000000000002a02 <+10754>:	addps  %xmm14,%xmm15
   0x0000000000002a42 <+10818>:	movaps 0x790(%r8),%xmm13
   0x0000000000002a4a <+10826>:	mulps  0x390(%rcx),%xmm13
   0x0000000000002a52 <+10834>:	addps  %xmm14,%xmm13
   0x0000000000002a92 <+10898>:	movaps 0x7d0(%r8),%xmm12
   0x0000000000002a9a <+10906>:	mulps  0x390(%rcx),%xmm12
   0x0000000000002aa2 <+10914>:	addps  %xmm14,%xmm12

853
854	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
855	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
856	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000029b6 <+10678>:	movaps 0x720(%r8),%xmm14
   0x00000000000029be <+10686>:	mulps  0x3a0(%rcx),%xmm14
   0x00000000000029c6 <+10694>:	addps  %xmm11,%xmm14
   0x0000000000002a06 <+10758>:	movaps 0x760(%r8),%xmm14
   0x0000000000002a0e <+10766>:	mulps  0x3a0(%rcx),%xmm14
   0x0000000000002a16 <+10774>:	addps  %xmm15,%xmm14
   0x0000000000002a56 <+10838>:	movaps 0x7a0(%r8),%xmm14
   0x0000000000002a5e <+10846>:	mulps  0x3a0(%rcx),%xmm14
   0x0000000000002a66 <+10854>:	addps  %xmm13,%xmm14
   0x0000000000002aa6 <+10918>:	movaps 0x7e0(%r8),%xmm14
   0x0000000000002aae <+10926>:	mulps  0x3a0(%rcx),%xmm14
   0x0000000000002ab6 <+10934>:	addps  %xmm12,%xmm14

857
858	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
859	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
860	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000029ca <+10698>:	movaps 0x730(%r8),%xmm11
   0x00000000000029d2 <+10706>:	mulps  0x3b0(%rcx),%xmm11
   0x00000000000029da <+10714>:	addps  %xmm14,%xmm11
   0x0000000000002a1a <+10778>:	movaps 0x770(%r8),%xmm15
   0x0000000000002a22 <+10786>:	mulps  0x3b0(%rcx),%xmm15
   0x0000000000002a2a <+10794>:	addps  %xmm14,%xmm15
   0x0000000000002a6a <+10858>:	movaps 0x7b0(%r8),%xmm13
   0x0000000000002a72 <+10866>:	mulps  0x3b0(%rcx),%xmm13
   0x0000000000002a7a <+10874>:	addps  %xmm14,%xmm13
   0x0000000000002aba <+10938>:	movaps 0x7f0(%r8),%xmm12
   0x0000000000002ac2 <+10946>:	mulps  0x3b0(%rcx),%xmm12
   0x0000000000002aca <+10954>:	addps  %xmm14,%xmm12

861	      }
862
863	      /* Store C(2,3) block. */
864	      for (i = 0; i < 4; i++)
865	      {
866	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000002ace <+10958>:	movaps 0x20(%rsp),%xmm14
   0x0000000000002ad4 <+10964>:	mulps  %xmm14,%xmm11
   0x0000000000002ae8 <+10984>:	mulps  %xmm14,%xmm15
   0x0000000000002afc <+11004>:	mulps  %xmm14,%xmm13
   0x0000000000002b10 <+11024>:	mulps  %xmm14,%xmm12

867	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_23]), C_row[i]);
   0x0000000000002ad8 <+10968>:	addps  0x180(%r9),%xmm11
   0x0000000000002aec <+10988>:	addps  0x190(%r9),%xmm15
   0x0000000000002b00 <+11008>:	addps  0x1a0(%r9),%xmm13
   0x0000000000002b14 <+11028>:	addps  0x1b0(%r9),%xmm12

868	        _mm_store_ps(&C[i*4+C_OFFSET_23], C_row[i]);
   0x0000000000002ae0 <+10976>:	movaps %xmm11,0x180(%r9)
   0x0000000000002af4 <+10996>:	movaps %xmm15,0x190(%r9)
   0x0000000000002b08 <+11016>:	movaps %xmm13,0x1a0(%r9)
   0x0000000000002b1c <+11036>:	movaps %xmm12,0x1b0(%r9)

869	      }
870	    }
871
872	    /* Reset C(2,4) matrix accumulators */
873	    C_row[0] = _mm_setzero_ps();
874	    C_row[1] = _mm_setzero_ps();
875	    C_row[2] = _mm_setzero_ps();
876	    C_row[3] = _mm_setzero_ps();
877
878	    if (norm_product[4][3] &&
   0x0000000000002b24 <+11044>:	comiss %xmm1,%xmm10
   0x0000000000002b28 <+11048>:	jb     0x30b0 <stream_kernel+12464>

879	        norm_product[5][7] &&
   0x0000000000002b2e <+11054>:	movss  0xf0(%rsp),%xmm10
   0x0000000000002b38 <+11064>:	comiss %xmm1,%xmm10
   0x0000000000002b3c <+11068>:	jb     0x30b0 <stream_kernel+12464>

880	        norm_product[6][11] &&
   0x0000000000002b42 <+11074>:	movss  0xd8(%rsp),%xmm10
   0x0000000000002b4c <+11084>:	comiss %xmm1,%xmm10
   0x0000000000002b50 <+11088>:	jb     0x30b0 <stream_kernel+12464>

881	        norm_product[7][15])
   0x0000000000002b56 <+11094>:	movss  0xd0(%rsp),%xmm10
   0x0000000000002b60 <+11104>:	comiss %xmm1,%xmm10
   0x0000000000002b64 <+11108>:	jb     0x30b0 <stream_kernel+12464>

882	    {
883	      /* A(2,1)*B(1,4) = C(2,4). */
884	      for (i = 0; i < 4; i++)
885	      {
886	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
887	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
888	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b6a <+11114>:	movaps 0x400(%r8),%xmm11
   0x0000000000002b8a <+11146>:	movaps 0x440(%r8),%xmm10
   0x0000000000002b9a <+11162>:	mulps  0xc0(%rcx),%xmm11
   0x0000000000002bba <+11194>:	mulps  0xc0(%rcx),%xmm10
   0x0000000000002bf6 <+11254>:	movaps 0x480(%r8),%xmm13
   0x0000000000002bfe <+11262>:	mulps  0xc0(%rcx),%xmm13
   0x0000000000002c32 <+11314>:	movaps 0x4c0(%r8),%xmm13
   0x0000000000002c3a <+11322>:	mulps  0xc0(%rcx),%xmm13

889
890	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
891	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
892	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b72 <+11122>:	movaps 0x410(%r8),%xmm15
   0x0000000000002b92 <+11154>:	movaps 0x450(%r8),%xmm12
   0x0000000000002ba2 <+11170>:	mulps  0xd0(%rcx),%xmm15
   0x0000000000002bc2 <+11202>:	mulps  0xd0(%rcx),%xmm12
   0x0000000000002bca <+11210>:	addps  %xmm11,%xmm15
   0x0000000000002bce <+11214>:	movaps 0x490(%r8),%xmm11
   0x0000000000002bd6 <+11222>:	mulps  0xd0(%rcx),%xmm11
   0x0000000000002c06 <+11270>:	addps  %xmm10,%xmm12
   0x0000000000002c2e <+11310>:	addps  %xmm13,%xmm11
   0x0000000000002c6e <+11374>:	movaps 0x4d0(%r8),%xmm15
   0x0000000000002c76 <+11382>:	mulps  0xd0(%rcx),%xmm15
   0x0000000000002c7e <+11390>:	addps  %xmm13,%xmm15

893
894	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
895	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
896	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b7a <+11130>:	movaps 0x420(%r8),%xmm13
   0x0000000000002baa <+11178>:	mulps  0xe0(%rcx),%xmm13
   0x0000000000002bde <+11230>:	addps  %xmm15,%xmm13
   0x0000000000002be2 <+11234>:	movaps 0x460(%r8),%xmm15
   0x0000000000002bea <+11242>:	mulps  0xe0(%rcx),%xmm15
   0x0000000000002c1a <+11290>:	addps  %xmm12,%xmm15
   0x0000000000002c46 <+11334>:	movaps 0x4a0(%r8),%xmm15
   0x0000000000002c4e <+11342>:	mulps  0xe0(%rcx),%xmm15
   0x0000000000002c56 <+11350>:	addps  %xmm11,%xmm15
   0x0000000000002c5a <+11354>:	movaps 0x4e0(%r8),%xmm11
   0x0000000000002c62 <+11362>:	mulps  0xe0(%rcx),%xmm11
   0x0000000000002c92 <+11410>:	addps  %xmm15,%xmm11

897
898	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
899	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
900	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b82 <+11138>:	movaps 0x430(%r8),%xmm14
   0x0000000000002bb2 <+11186>:	mulps  0xf0(%rcx),%xmm14
   0x0000000000002bf2 <+11250>:	addps  %xmm13,%xmm14
   0x0000000000002c0a <+11274>:	movaps 0x4b0(%r8),%xmm10
   0x0000000000002c12 <+11282>:	mulps  0xf0(%rcx),%xmm10
   0x0000000000002c1e <+11294>:	movaps 0x470(%r8),%xmm12
   0x0000000000002c26 <+11302>:	mulps  0xf0(%rcx),%xmm12
   0x0000000000002c42 <+11330>:	addps  %xmm15,%xmm12
   0x0000000000002c6a <+11370>:	addps  %xmm15,%xmm10
   0x0000000000002c82 <+11394>:	movaps 0x4f0(%r8),%xmm13
   0x0000000000002c8a <+11402>:	mulps  0xf0(%rcx),%xmm13
   0x0000000000002ca6 <+11430>:	addps  %xmm11,%xmm13

901	      }
902
903	      /* A(2,2)*B(2,4) = C(2,4). */
904	      for (i = 0; i < 4; i++)
905	      {
906	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
907	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
908	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002c96 <+11414>:	movaps 0x500(%r8),%xmm15
   0x0000000000002c9e <+11422>:	mulps  0x1c0(%rcx),%xmm15
   0x0000000000002cba <+11450>:	addps  %xmm14,%xmm15
   0x0000000000002cd2 <+11474>:	movaps 0x540(%r8),%xmm15
   0x0000000000002cda <+11482>:	mulps  0x1c0(%rcx),%xmm15
   0x0000000000002cf6 <+11510>:	addps  %xmm12,%xmm15
   0x0000000000002d36 <+11574>:	movaps 0x580(%r8),%xmm12
   0x0000000000002d3e <+11582>:	mulps  0x1c0(%rcx),%xmm12
   0x0000000000002d5a <+11610>:	addps  %xmm10,%xmm12
   0x0000000000002d86 <+11654>:	movaps 0x5c0(%r8),%xmm10
   0x0000000000002d8e <+11662>:	mulps  0x1c0(%rcx),%xmm10
   0x0000000000002daa <+11690>:	addps  %xmm13,%xmm10

909
910	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
911	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
912	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002cbe <+11454>:	movaps 0x510(%r8),%xmm14
   0x0000000000002cc6 <+11462>:	mulps  0x1d0(%rcx),%xmm14
   0x0000000000002cce <+11470>:	addps  %xmm15,%xmm14
   0x0000000000002cfa <+11514>:	movaps 0x550(%r8),%xmm12
   0x0000000000002d02 <+11522>:	mulps  0x1d0(%rcx),%xmm12
   0x0000000000002d1e <+11550>:	addps  %xmm15,%xmm12
   0x0000000000002d5e <+11614>:	movaps 0x590(%r8),%xmm10
   0x0000000000002d66 <+11622>:	mulps  0x1d0(%rcx),%xmm10
   0x0000000000002d6e <+11630>:	addps  %xmm12,%xmm10
   0x0000000000002dae <+11694>:	movaps 0x5d0(%r8),%xmm13
   0x0000000000002db6 <+11702>:	mulps  0x1d0(%rcx),%xmm13
   0x0000000000002dbe <+11710>:	addps  %xmm10,%xmm13

913
914	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
915	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
916	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002caa <+11434>:	movaps 0x520(%r8),%xmm11
   0x0000000000002cb2 <+11442>:	mulps  0x1e0(%rcx),%xmm11
   0x0000000000002ce2 <+11490>:	addps  %xmm14,%xmm11
   0x0000000000002d0e <+11534>:	movaps 0x560(%r8),%xmm11
   0x0000000000002d16 <+11542>:	mulps  0x1e0(%rcx),%xmm11
   0x0000000000002d32 <+11570>:	addps  %xmm12,%xmm11
   0x0000000000002d4a <+11594>:	movaps 0x5a0(%r8),%xmm11
   0x0000000000002d52 <+11602>:	mulps  0x1e0(%rcx),%xmm11
   0x0000000000002d82 <+11650>:	addps  %xmm10,%xmm11
   0x0000000000002d9a <+11674>:	movaps 0x5e0(%r8),%xmm11
   0x0000000000002da2 <+11682>:	mulps  0x1e0(%rcx),%xmm11
   0x0000000000002dd2 <+11730>:	addps  %xmm13,%xmm11

917
918	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
919	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
920	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002ce6 <+11494>:	movaps 0x530(%r8),%xmm14
   0x0000000000002cee <+11502>:	mulps  0x1f0(%rcx),%xmm14
   0x0000000000002d0a <+11530>:	addps  %xmm11,%xmm14
   0x0000000000002d22 <+11554>:	movaps 0x570(%r8),%xmm15
   0x0000000000002d2a <+11562>:	mulps  0x1f0(%rcx),%xmm15
   0x0000000000002d46 <+11590>:	addps  %xmm11,%xmm15
   0x0000000000002d72 <+11634>:	movaps 0x5b0(%r8),%xmm12
   0x0000000000002d7a <+11642>:	mulps  0x1f0(%rcx),%xmm12
   0x0000000000002d96 <+11670>:	addps  %xmm11,%xmm12
   0x0000000000002dc2 <+11714>:	movaps 0x5f0(%r8),%xmm10
   0x0000000000002dca <+11722>:	mulps  0x1f0(%rcx),%xmm10
   0x0000000000002de6 <+11750>:	addps  %xmm11,%xmm10

921	      }
922
923	      /* A(2,3)*B(3,4) = C(2,4). */
924	      for (i = 0; i < 4; i++)
925	      {
926	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
927	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
928	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002dd6 <+11734>:	movaps 0x600(%r8),%xmm13
   0x0000000000002dde <+11742>:	mulps  0x2c0(%rcx),%xmm13
   0x0000000000002dfa <+11770>:	addps  %xmm14,%xmm13
   0x0000000000002e26 <+11814>:	movaps 0x640(%r8),%xmm14
   0x0000000000002e2e <+11822>:	mulps  0x2c0(%rcx),%xmm14
   0x0000000000002e4a <+11850>:	addps  %xmm15,%xmm14
   0x0000000000002e76 <+11894>:	movaps 0x680(%r8),%xmm15
   0x0000000000002e7e <+11902>:	mulps  0x2c0(%rcx),%xmm15
   0x0000000000002e9a <+11930>:	addps  %xmm12,%xmm15
   0x0000000000002eb2 <+11954>:	movaps 0x6c0(%r8),%xmm15
   0x0000000000002eba <+11962>:	mulps  0x2c0(%rcx),%xmm15
   0x0000000000002ed6 <+11990>:	addps  %xmm10,%xmm15

929
930	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
931	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
932	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002dfe <+11774>:	movaps 0x610(%r8),%xmm14
   0x0000000000002e06 <+11782>:	mulps  0x2d0(%rcx),%xmm14
   0x0000000000002e0e <+11790>:	addps  %xmm13,%xmm14
   0x0000000000002e4e <+11854>:	movaps 0x650(%r8),%xmm15
   0x0000000000002e56 <+11862>:	mulps  0x2d0(%rcx),%xmm15
   0x0000000000002e5e <+11870>:	addps  %xmm14,%xmm15
   0x0000000000002e9e <+11934>:	movaps 0x690(%r8),%xmm12
   0x0000000000002ea6 <+11942>:	mulps  0x2d0(%rcx),%xmm12
   0x0000000000002eae <+11950>:	addps  %xmm15,%xmm12
   0x0000000000002eda <+11994>:	movaps 0x6d0(%r8),%xmm10
   0x0000000000002ee2 <+12002>:	mulps  0x2d0(%rcx),%xmm10
   0x0000000000002efe <+12030>:	addps  %xmm15,%xmm10

933
934	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
935	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
936	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002e12 <+11794>:	movaps 0x620(%r8),%xmm13
   0x0000000000002e1a <+11802>:	mulps  0x2e0(%rcx),%xmm13
   0x0000000000002e22 <+11810>:	addps  %xmm14,%xmm13
   0x0000000000002e3a <+11834>:	movaps 0x660(%r8),%xmm13
   0x0000000000002e42 <+11842>:	mulps  0x2e0(%rcx),%xmm13
   0x0000000000002e72 <+11890>:	addps  %xmm15,%xmm13
   0x0000000000002e8a <+11914>:	movaps 0x6a0(%r8),%xmm13
   0x0000000000002e92 <+11922>:	mulps  0x2e0(%rcx),%xmm13
   0x0000000000002ec2 <+11970>:	addps  %xmm12,%xmm13
   0x0000000000002eee <+12014>:	movaps 0x6e0(%r8),%xmm13
   0x0000000000002ef6 <+12022>:	mulps  0x2e0(%rcx),%xmm13
   0x0000000000002f12 <+12050>:	addps  %xmm10,%xmm13

937
938	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
939	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
940	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002dea <+11754>:	movaps 0x630(%r8),%xmm11
   0x0000000000002df2 <+11762>:	mulps  0x2f0(%rcx),%xmm11
   0x0000000000002e36 <+11830>:	addps  %xmm13,%xmm11
   0x0000000000002e62 <+11874>:	movaps 0x670(%r8),%xmm14
   0x0000000000002e6a <+11882>:	mulps  0x2f0(%rcx),%xmm14
   0x0000000000002e86 <+11910>:	addps  %xmm13,%xmm14
   0x0000000000002ec6 <+11974>:	movaps 0x6b0(%r8),%xmm12
   0x0000000000002ece <+11982>:	mulps  0x2f0(%rcx),%xmm12
   0x0000000000002eea <+12010>:	addps  %xmm13,%xmm12
   0x0000000000002f16 <+12054>:	movaps 0x6f0(%r8),%xmm10
   0x0000000000002f1e <+12062>:	mulps  0x2f0(%rcx),%xmm10
   0x0000000000002f3a <+12090>:	addps  %xmm13,%xmm10

941	      }
942
943	      /* A(2,4)*B(4,4) = C(2,4). */
944	      for (i = 0; i < 4; i++)
945	      {
946	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
947	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
948	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002f02 <+12034>:	movaps 0x700(%r8),%xmm15
   0x0000000000002f0a <+12042>:	mulps  0x3c0(%rcx),%xmm15
   0x0000000000002f26 <+12070>:	addps  %xmm11,%xmm15
   0x0000000000002f66 <+12134>:	movaps 0x740(%r8),%xmm11
   0x0000000000002f6e <+12142>:	mulps  0x3c0(%rcx),%xmm11
   0x0000000000002f8a <+12170>:	addps  %xmm14,%xmm11
   0x0000000000002fb6 <+12214>:	movaps 0x780(%r8),%xmm14
   0x0000000000002fbe <+12222>:	mulps  0x3c0(%rcx),%xmm14
   0x0000000000002fda <+12250>:	addps  %xmm12,%xmm14
   0x0000000000003006 <+12294>:	movaps 0x7c0(%r8),%xmm12
   0x000000000000300e <+12302>:	mulps  0x3c0(%rcx),%xmm12
   0x000000000000302a <+12330>:	addps  %xmm10,%xmm12

949
950	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
951	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
952	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002f2a <+12074>:	movaps 0x710(%r8),%xmm11
   0x0000000000002f32 <+12082>:	mulps  0x3d0(%rcx),%xmm11
   0x0000000000002f4e <+12110>:	addps  %xmm15,%xmm11
   0x0000000000002f8e <+12174>:	movaps 0x750(%r8),%xmm14
   0x0000000000002f96 <+12182>:	mulps  0x3d0(%rcx),%xmm14
   0x0000000000002f9e <+12190>:	addps  %xmm11,%xmm14
   0x0000000000002fde <+12254>:	movaps 0x790(%r8),%xmm12
   0x0000000000002fe6 <+12262>:	mulps  0x3d0(%rcx),%xmm12
   0x0000000000002fee <+12270>:	addps  %xmm14,%xmm12
   0x000000000000302e <+12334>:	movaps 0x7d0(%r8),%xmm10
   0x0000000000003036 <+12342>:	mulps  0x3d0(%rcx),%xmm10
   0x000000000000303e <+12350>:	addps  %xmm12,%xmm10

953
954	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
955	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
956	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002f3e <+12094>:	movaps 0x720(%r8),%xmm13
   0x0000000000002f46 <+12102>:	mulps  0x3e0(%rcx),%xmm13
   0x0000000000002f62 <+12130>:	addps  %xmm11,%xmm13
   0x0000000000002fa2 <+12194>:	movaps 0x760(%r8),%xmm11
   0x0000000000002faa <+12202>:	mulps  0x3e0(%rcx),%xmm11
   0x0000000000002fb2 <+12210>:	addps  %xmm14,%xmm11
   0x0000000000002ff2 <+12274>:	movaps 0x7a0(%r8),%xmm14
   0x0000000000002ffa <+12282>:	mulps  0x3e0(%rcx),%xmm14
   0x0000000000003002 <+12290>:	addps  %xmm12,%xmm14
   0x0000000000003042 <+12354>:	movaps 0x7e0(%r8),%xmm12
   0x000000000000304a <+12362>:	mulps  0x3e0(%rcx),%xmm12
   0x0000000000003052 <+12370>:	addps  %xmm10,%xmm12

957
958	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
959	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
960	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002f52 <+12114>:	movaps 0x730(%r8),%xmm15
   0x0000000000002f5a <+12122>:	mulps  0x3f0(%rcx),%xmm15
   0x0000000000002f76 <+12150>:	addps  %xmm13,%xmm15
   0x0000000000002f7a <+12154>:	movaps 0x770(%r8),%xmm13
   0x0000000000002f82 <+12162>:	mulps  0x3f0(%rcx),%xmm13
   0x0000000000002fc6 <+12230>:	addps  %xmm11,%xmm13
   0x0000000000002fca <+12234>:	movaps 0x7b0(%r8),%xmm11
   0x0000000000002fd2 <+12242>:	mulps  0x3f0(%rcx),%xmm11
   0x0000000000003016 <+12310>:	addps  %xmm14,%xmm11
   0x000000000000301a <+12314>:	movaps 0x7f0(%r8),%xmm14
   0x0000000000003022 <+12322>:	mulps  0x3f0(%rcx),%xmm14
   0x0000000000003098 <+12440>:	addps  %xmm12,%xmm14

961	      }
962
963	      /* Store C(2,4) block. */
964	      for (i = 0; i < 4; i++)
965	      {
966	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000003056 <+12374>:	movaps 0x20(%rsp),%xmm10
   0x000000000000305c <+12380>:	mulps  %xmm10,%xmm15
   0x0000000000003070 <+12400>:	mulps  %xmm10,%xmm13
   0x0000000000003084 <+12420>:	mulps  %xmm10,%xmm11
   0x000000000000309c <+12444>:	mulps  %xmm10,%xmm14

967	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_24]), C_row[i]);
   0x0000000000003060 <+12384>:	addps  0x1c0(%r9),%xmm15
   0x0000000000003074 <+12404>:	addps  0x1d0(%r9),%xmm13
   0x0000000000003088 <+12424>:	addps  0x1e0(%r9),%xmm11
   0x00000000000030a0 <+12448>:	addps  0x1f0(%r9),%xmm14

968	        _mm_store_ps(&C[i*4+C_OFFSET_24], C_row[i]);
   0x0000000000003068 <+12392>:	movaps %xmm15,0x1c0(%r9)
   0x000000000000307c <+12412>:	movaps %xmm13,0x1d0(%r9)
   0x0000000000003090 <+12432>:	movaps %xmm11,0x1e0(%r9)
   0x00000000000030a8 <+12456>:	movaps %xmm14,0x1f0(%r9)

969	      }
970	    }
971
972	    /* Reset C(3,1) matrix accumulators */
973	    C_row[0] = _mm_setzero_ps();
974	    C_row[1] = _mm_setzero_ps();
975	    C_row[2] = _mm_setzero_ps();
976	    C_row[3] = _mm_setzero_ps();
977
978	    if (norm_product[8][0] &&
   0x00000000000030b0 <+12464>:	movss  0x160(%rsp),%xmm10
   0x00000000000030ba <+12474>:	comiss %xmm1,%xmm10
   0x00000000000030be <+12478>:	jb     0x3612 <stream_kernel+13842>

979	        norm_product[9][4] &&
   0x00000000000030c4 <+12484>:	movss  0xc8(%rsp),%xmm10
   0x00000000000030ce <+12494>:	comiss %xmm1,%xmm10
   0x00000000000030d2 <+12498>:	jb     0x3612 <stream_kernel+13842>

980	        norm_product[10][8] &&
   0x00000000000030d8 <+12504>:	movss  0xc0(%rsp),%xmm10
   0x00000000000030e2 <+12514>:	comiss %xmm1,%xmm10
   0x00000000000030e6 <+12518>:	jb     0x3612 <stream_kernel+13842>

981	        norm_product[11][12])
   0x00000000000030ec <+12524>:	movss  0xa8(%rsp),%xmm10
   0x00000000000030f6 <+12534>:	comiss %xmm1,%xmm10
   0x00000000000030fa <+12538>:	jb     0x3612 <stream_kernel+13842>

982	    {
983	      /* A(3,1)*B(1,1) = C(3,1). */
984	      for (i = 0; i < 4; i++)
985	      {
986	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
987	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
988	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003100 <+12544>:	movaps 0x800(%r8),%xmm11
   0x0000000000003120 <+12576>:	movaps 0x840(%r8),%xmm10
   0x0000000000003130 <+12592>:	mulps  (%rcx),%xmm11
   0x0000000000003143 <+12611>:	mulps  (%rcx),%xmm10
   0x0000000000003172 <+12658>:	movaps 0x880(%r8),%xmm13
   0x000000000000317a <+12666>:	mulps  (%rcx),%xmm13
   0x00000000000031a4 <+12708>:	movaps 0x8c0(%r8),%xmm13
   0x00000000000031ac <+12716>:	mulps  (%rcx),%xmm13

989
990	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
991	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
992	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003108 <+12552>:	movaps 0x810(%r8),%xmm15
   0x0000000000003128 <+12584>:	movaps 0x850(%r8),%xmm12
   0x0000000000003134 <+12596>:	mulps  0x10(%rcx),%xmm15
   0x0000000000003147 <+12615>:	mulps  0x10(%rcx),%xmm12
   0x000000000000314c <+12620>:	addps  %xmm11,%xmm15
   0x0000000000003150 <+12624>:	movaps 0x890(%r8),%xmm11
   0x0000000000003158 <+12632>:	mulps  0x10(%rcx),%xmm11
   0x000000000000317e <+12670>:	addps  %xmm10,%xmm12
   0x00000000000031a0 <+12704>:	addps  %xmm13,%xmm11
   0x00000000000031d6 <+12758>:	movaps 0x8d0(%r8),%xmm15
   0x00000000000031de <+12766>:	mulps  0x10(%rcx),%xmm15
   0x00000000000031e3 <+12771>:	addps  %xmm13,%xmm15

993
994	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
995	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
996	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003110 <+12560>:	movaps 0x820(%r8),%xmm13
   0x0000000000003139 <+12601>:	mulps  0x20(%rcx),%xmm13
   0x000000000000315d <+12637>:	addps  %xmm15,%xmm13
   0x0000000000003161 <+12641>:	movaps 0x860(%r8),%xmm15
   0x0000000000003169 <+12649>:	mulps  0x20(%rcx),%xmm15
   0x000000000000318f <+12687>:	addps  %xmm12,%xmm15
   0x00000000000031b4 <+12724>:	movaps 0x8a0(%r8),%xmm15
   0x00000000000031bc <+12732>:	mulps  0x20(%rcx),%xmm15
   0x00000000000031c1 <+12737>:	addps  %xmm11,%xmm15
   0x00000000000031c5 <+12741>:	movaps 0x8e0(%r8),%xmm11
   0x00000000000031cd <+12749>:	mulps  0x20(%rcx),%xmm11
   0x00000000000031f4 <+12788>:	addps  %xmm15,%xmm11

997
998	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
999	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
1000	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003118 <+12568>:	movaps 0x830(%r8),%xmm14
   0x000000000000313e <+12606>:	mulps  0x30(%rcx),%xmm14
   0x000000000000316e <+12654>:	addps  %xmm13,%xmm14
   0x0000000000003182 <+12674>:	movaps 0x8b0(%r8),%xmm10
   0x000000000000318a <+12682>:	mulps  0x30(%rcx),%xmm10
   0x0000000000003193 <+12691>:	movaps 0x870(%r8),%xmm12
   0x000000000000319b <+12699>:	mulps  0x30(%rcx),%xmm12
   0x00000000000031b0 <+12720>:	addps  %xmm15,%xmm12
   0x00000000000031d2 <+12754>:	addps  %xmm15,%xmm10
   0x00000000000031e7 <+12775>:	movaps 0x8f0(%r8),%xmm13
   0x00000000000031ef <+12783>:	mulps  0x30(%rcx),%xmm13
   0x0000000000003208 <+12808>:	addps  %xmm11,%xmm13

1001	      }
1002
1003	      /* A(3,2)*B(2,1) = C(3,1). */
1004	      for (i = 0; i < 4; i++)
1005	      {
1006	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1007	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
1008	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000031f8 <+12792>:	movaps 0x900(%r8),%xmm15
   0x0000000000003200 <+12800>:	mulps  0x100(%rcx),%xmm15
   0x000000000000321c <+12828>:	addps  %xmm14,%xmm15
   0x0000000000003234 <+12852>:	movaps 0x940(%r8),%xmm15
   0x000000000000323c <+12860>:	mulps  0x100(%rcx),%xmm15
   0x0000000000003258 <+12888>:	addps  %xmm12,%xmm15
   0x0000000000003298 <+12952>:	movaps 0x980(%r8),%xmm12
   0x00000000000032a0 <+12960>:	mulps  0x100(%rcx),%xmm12
   0x00000000000032bc <+12988>:	addps  %xmm10,%xmm12
   0x00000000000032e8 <+13032>:	movaps 0x9c0(%r8),%xmm10
   0x00000000000032f0 <+13040>:	mulps  0x100(%rcx),%xmm10
   0x000000000000330c <+13068>:	addps  %xmm13,%xmm10

1009
1010	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1011	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
1012	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003220 <+12832>:	movaps 0x910(%r8),%xmm14
   0x0000000000003228 <+12840>:	mulps  0x110(%rcx),%xmm14
   0x0000000000003230 <+12848>:	addps  %xmm15,%xmm14
   0x000000000000325c <+12892>:	movaps 0x950(%r8),%xmm12
   0x0000000000003264 <+12900>:	mulps  0x110(%rcx),%xmm12
   0x0000000000003280 <+12928>:	addps  %xmm15,%xmm12
   0x00000000000032c0 <+12992>:	movaps 0x990(%r8),%xmm10
   0x00000000000032c8 <+13000>:	mulps  0x110(%rcx),%xmm10
   0x00000000000032d0 <+13008>:	addps  %xmm12,%xmm10
   0x0000000000003310 <+13072>:	movaps 0x9d0(%r8),%xmm13
   0x0000000000003318 <+13080>:	mulps  0x110(%rcx),%xmm13
   0x0000000000003320 <+13088>:	addps  %xmm10,%xmm13

1013
1014	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1015	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
1016	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000320c <+12812>:	movaps 0x920(%r8),%xmm11
   0x0000000000003214 <+12820>:	mulps  0x120(%rcx),%xmm11
   0x0000000000003244 <+12868>:	addps  %xmm14,%xmm11
   0x0000000000003270 <+12912>:	movaps 0x960(%r8),%xmm11
   0x0000000000003278 <+12920>:	mulps  0x120(%rcx),%xmm11
   0x0000000000003294 <+12948>:	addps  %xmm12,%xmm11
   0x00000000000032ac <+12972>:	movaps 0x9a0(%r8),%xmm11
   0x00000000000032b4 <+12980>:	mulps  0x120(%rcx),%xmm11
   0x00000000000032e4 <+13028>:	addps  %xmm10,%xmm11
   0x00000000000032fc <+13052>:	movaps 0x9e0(%r8),%xmm11
   0x0000000000003304 <+13060>:	mulps  0x120(%rcx),%xmm11
   0x0000000000003334 <+13108>:	addps  %xmm13,%xmm11

1017
1018	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1019	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
1020	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003248 <+12872>:	movaps 0x930(%r8),%xmm14
   0x0000000000003250 <+12880>:	mulps  0x130(%rcx),%xmm14
   0x000000000000326c <+12908>:	addps  %xmm11,%xmm14
   0x0000000000003284 <+12932>:	movaps 0x970(%r8),%xmm15
   0x000000000000328c <+12940>:	mulps  0x130(%rcx),%xmm15
   0x00000000000032a8 <+12968>:	addps  %xmm11,%xmm15
   0x00000000000032d4 <+13012>:	movaps 0x9b0(%r8),%xmm12
   0x00000000000032dc <+13020>:	mulps  0x130(%rcx),%xmm12
   0x00000000000032f8 <+13048>:	addps  %xmm11,%xmm12
   0x0000000000003324 <+13092>:	movaps 0x9f0(%r8),%xmm10
   0x000000000000332c <+13100>:	mulps  0x130(%rcx),%xmm10
   0x0000000000003348 <+13128>:	addps  %xmm11,%xmm10

1021	      }
1022
1023	      /* A(3,3)*B(3,1) = C(3,1). */
1024	      for (i = 0; i < 4; i++)
1025	      {
1026	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1027	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
1028	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003338 <+13112>:	movaps 0xa00(%r8),%xmm13
   0x0000000000003340 <+13120>:	mulps  0x200(%rcx),%xmm13
   0x000000000000335c <+13148>:	addps  %xmm14,%xmm13
   0x0000000000003388 <+13192>:	movaps 0xa40(%r8),%xmm14
   0x0000000000003390 <+13200>:	mulps  0x200(%rcx),%xmm14
   0x00000000000033ac <+13228>:	addps  %xmm15,%xmm14
   0x00000000000033d8 <+13272>:	movaps 0xa80(%r8),%xmm15
   0x00000000000033e0 <+13280>:	mulps  0x200(%rcx),%xmm15
   0x00000000000033fc <+13308>:	addps  %xmm12,%xmm15
   0x0000000000003414 <+13332>:	movaps 0xac0(%r8),%xmm15
   0x000000000000341c <+13340>:	mulps  0x200(%rcx),%xmm15
   0x0000000000003438 <+13368>:	addps  %xmm10,%xmm15

1029
1030	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1031	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
1032	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003360 <+13152>:	movaps 0xa10(%r8),%xmm14
   0x0000000000003368 <+13160>:	mulps  0x210(%rcx),%xmm14
   0x0000000000003370 <+13168>:	addps  %xmm13,%xmm14
   0x00000000000033b0 <+13232>:	movaps 0xa50(%r8),%xmm15
   0x00000000000033b8 <+13240>:	mulps  0x210(%rcx),%xmm15
   0x00000000000033c0 <+13248>:	addps  %xmm14,%xmm15
   0x0000000000003400 <+13312>:	movaps 0xa90(%r8),%xmm12
   0x0000000000003408 <+13320>:	mulps  0x210(%rcx),%xmm12
   0x0000000000003410 <+13328>:	addps  %xmm15,%xmm12
   0x000000000000343c <+13372>:	movaps 0xad0(%r8),%xmm10
   0x0000000000003444 <+13380>:	mulps  0x210(%rcx),%xmm10
   0x0000000000003460 <+13408>:	addps  %xmm15,%xmm10

1033
1034	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1035	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
1036	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003374 <+13172>:	movaps 0xa20(%r8),%xmm13
   0x000000000000337c <+13180>:	mulps  0x220(%rcx),%xmm13
   0x0000000000003384 <+13188>:	addps  %xmm14,%xmm13
   0x000000000000339c <+13212>:	movaps 0xa60(%r8),%xmm13
   0x00000000000033a4 <+13220>:	mulps  0x220(%rcx),%xmm13
   0x00000000000033d4 <+13268>:	addps  %xmm15,%xmm13
   0x00000000000033ec <+13292>:	movaps 0xaa0(%r8),%xmm13
   0x00000000000033f4 <+13300>:	mulps  0x220(%rcx),%xmm13
   0x0000000000003424 <+13348>:	addps  %xmm12,%xmm13
   0x0000000000003450 <+13392>:	movaps 0xae0(%r8),%xmm13
   0x0000000000003458 <+13400>:	mulps  0x220(%rcx),%xmm13
   0x0000000000003474 <+13428>:	addps  %xmm10,%xmm13

1037
1038	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1039	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
1040	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000334c <+13132>:	movaps 0xa30(%r8),%xmm11
   0x0000000000003354 <+13140>:	mulps  0x230(%rcx),%xmm11
   0x0000000000003398 <+13208>:	addps  %xmm13,%xmm11
   0x00000000000033c4 <+13252>:	movaps 0xa70(%r8),%xmm14
   0x00000000000033cc <+13260>:	mulps  0x230(%rcx),%xmm14
   0x00000000000033e8 <+13288>:	addps  %xmm13,%xmm14
   0x0000000000003428 <+13352>:	movaps 0xab0(%r8),%xmm12
   0x0000000000003430 <+13360>:	mulps  0x230(%rcx),%xmm12
   0x000000000000344c <+13388>:	addps  %xmm13,%xmm12
   0x0000000000003478 <+13432>:	movaps 0xaf0(%r8),%xmm10
   0x0000000000003480 <+13440>:	mulps  0x230(%rcx),%xmm10
   0x000000000000349c <+13468>:	addps  %xmm13,%xmm10

1041	      }
1042
1043	      /* A(3,4)*B(4,1) = C(3,1). */
1044	      for (i = 0; i < 4; i++)
1045	      {
1046	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1047	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
1048	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003464 <+13412>:	movaps 0xb00(%r8),%xmm15
   0x000000000000346c <+13420>:	mulps  0x300(%rcx),%xmm15
   0x0000000000003488 <+13448>:	addps  %xmm11,%xmm15
   0x00000000000034c8 <+13512>:	movaps 0xb40(%r8),%xmm11
   0x00000000000034d0 <+13520>:	mulps  0x300(%rcx),%xmm11
   0x00000000000034ec <+13548>:	addps  %xmm14,%xmm11
   0x0000000000003518 <+13592>:	movaps 0xb80(%r8),%xmm14
   0x0000000000003520 <+13600>:	mulps  0x300(%rcx),%xmm14
   0x000000000000353c <+13628>:	addps  %xmm12,%xmm14
   0x0000000000003568 <+13672>:	movaps 0xbc0(%r8),%xmm12
   0x0000000000003570 <+13680>:	mulps  0x300(%rcx),%xmm12
   0x000000000000358c <+13708>:	addps  %xmm10,%xmm12

1049
1050	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1051	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
1052	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000348c <+13452>:	movaps 0xb10(%r8),%xmm11
   0x0000000000003494 <+13460>:	mulps  0x310(%rcx),%xmm11
   0x00000000000034b0 <+13488>:	addps  %xmm15,%xmm11
   0x00000000000034f0 <+13552>:	movaps 0xb50(%r8),%xmm14
   0x00000000000034f8 <+13560>:	mulps  0x310(%rcx),%xmm14
   0x0000000000003500 <+13568>:	addps  %xmm11,%xmm14
   0x0000000000003540 <+13632>:	movaps 0xb90(%r8),%xmm12
   0x0000000000003548 <+13640>:	mulps  0x310(%rcx),%xmm12
   0x0000000000003550 <+13648>:	addps  %xmm14,%xmm12
   0x0000000000003590 <+13712>:	movaps 0xbd0(%r8),%xmm10
   0x0000000000003598 <+13720>:	mulps  0x310(%rcx),%xmm10
   0x00000000000035a0 <+13728>:	addps  %xmm12,%xmm10

1053
1054	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1055	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
1056	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000034a0 <+13472>:	movaps 0xb20(%r8),%xmm13
   0x00000000000034a8 <+13480>:	mulps  0x320(%rcx),%xmm13
   0x00000000000034c4 <+13508>:	addps  %xmm11,%xmm13
   0x0000000000003504 <+13572>:	movaps 0xb60(%r8),%xmm11
   0x000000000000350c <+13580>:	mulps  0x320(%rcx),%xmm11
   0x0000000000003514 <+13588>:	addps  %xmm14,%xmm11
   0x0000000000003554 <+13652>:	movaps 0xba0(%r8),%xmm14
   0x000000000000355c <+13660>:	mulps  0x320(%rcx),%xmm14
   0x0000000000003564 <+13668>:	addps  %xmm12,%xmm14
   0x00000000000035a4 <+13732>:	movaps 0xbe0(%r8),%xmm12
   0x00000000000035ac <+13740>:	mulps  0x320(%rcx),%xmm12
   0x00000000000035b4 <+13748>:	addps  %xmm10,%xmm12

1057
1058	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1059	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
1060	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000034b4 <+13492>:	movaps 0xb30(%r8),%xmm15
   0x00000000000034bc <+13500>:	mulps  0x330(%rcx),%xmm15
   0x00000000000034d8 <+13528>:	addps  %xmm13,%xmm15
   0x00000000000034dc <+13532>:	movaps 0xb70(%r8),%xmm13
   0x00000000000034e4 <+13540>:	mulps  0x330(%rcx),%xmm13
   0x0000000000003528 <+13608>:	addps  %xmm11,%xmm13
   0x000000000000352c <+13612>:	movaps 0xbb0(%r8),%xmm11
   0x0000000000003534 <+13620>:	mulps  0x330(%rcx),%xmm11
   0x0000000000003578 <+13688>:	addps  %xmm14,%xmm11
   0x000000000000357c <+13692>:	movaps 0xbf0(%r8),%xmm14
   0x0000000000003584 <+13700>:	mulps  0x330(%rcx),%xmm14
   0x00000000000035fa <+13818>:	addps  %xmm12,%xmm14

1061	      }
1062
1063	      /* Store C(3,1) block. */
1064	      for (i = 0; i < 4; i++)
1065	      {
1066	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000035b8 <+13752>:	movaps 0x20(%rsp),%xmm10
   0x00000000000035be <+13758>:	mulps  %xmm10,%xmm15
   0x00000000000035d2 <+13778>:	mulps  %xmm10,%xmm13
   0x00000000000035e6 <+13798>:	mulps  %xmm10,%xmm11
   0x00000000000035fe <+13822>:	mulps  %xmm10,%xmm14

1067	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_31]), C_row[i]);
   0x00000000000035c2 <+13762>:	addps  0x200(%r9),%xmm15
   0x00000000000035d6 <+13782>:	addps  0x210(%r9),%xmm13
   0x00000000000035ea <+13802>:	addps  0x220(%r9),%xmm11
   0x0000000000003602 <+13826>:	addps  0x230(%r9),%xmm14

1068	        _mm_store_ps(&C[i*4+C_OFFSET_31], C_row[i]);
   0x00000000000035ca <+13770>:	movaps %xmm15,0x200(%r9)
   0x00000000000035de <+13790>:	movaps %xmm13,0x210(%r9)
   0x00000000000035f2 <+13810>:	movaps %xmm11,0x220(%r9)
   0x000000000000360a <+13834>:	movaps %xmm14,0x230(%r9)

1069	      }
1070	    }
1071
1072	    /* Reset C(3,2) matrix accumulators */
1073	    C_row[0] = _mm_setzero_ps();
1074	    C_row[1] = _mm_setzero_ps();
1075	    C_row[2] = _mm_setzero_ps();
1076	    C_row[3] = _mm_setzero_ps();
1077
1078	    if (norm_product[8][1] &&
   0x0000000000003612 <+13842>:	movss  0x128(%rsp),%xmm10
   0x000000000000361c <+13852>:	comiss %xmm1,%xmm10
   0x0000000000003620 <+13856>:	jb     0x3b78 <stream_kernel+15224>

1079	        norm_product[9][5] &&
   0x0000000000003626 <+13862>:	movss  0xb0(%rsp),%xmm10
   0x0000000000003630 <+13872>:	comiss %xmm1,%xmm10
   0x0000000000003634 <+13876>:	jb     0x3b78 <stream_kernel+15224>

1080	        norm_product[10][9] &&
   0x000000000000363a <+13882>:	movss  0x98(%rsp),%xmm10
   0x0000000000003644 <+13892>:	comiss %xmm1,%xmm10
   0x0000000000003648 <+13896>:	jb     0x3b78 <stream_kernel+15224>

1081	        norm_product[11][13])
   0x000000000000364e <+13902>:	movss  0x88(%rsp),%xmm10
   0x0000000000003658 <+13912>:	comiss %xmm1,%xmm10
   0x000000000000365c <+13916>:	jb     0x3b78 <stream_kernel+15224>

1082	    {
1083	      /* A(3,1)*B(1,2) = C(3,2). */
1084	      for (i = 0; i < 4; i++)
1085	      {
1086	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
1087	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
1088	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003662 <+13922>:	movaps 0x800(%r8),%xmm11
   0x0000000000003682 <+13954>:	movaps 0x840(%r8),%xmm10
   0x0000000000003692 <+13970>:	mulps  0x40(%rcx),%xmm11
   0x00000000000036a6 <+13990>:	mulps  0x40(%rcx),%xmm10
   0x00000000000036d6 <+14038>:	movaps 0x880(%r8),%xmm13
   0x00000000000036de <+14046>:	mulps  0x40(%rcx),%xmm13
   0x0000000000003709 <+14089>:	movaps 0x8c0(%r8),%xmm13
   0x0000000000003711 <+14097>:	mulps  0x40(%rcx),%xmm13

1089
1090	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
1091	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
1092	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000366a <+13930>:	movaps 0x810(%r8),%xmm15
   0x000000000000368a <+13962>:	movaps 0x850(%r8),%xmm12
   0x0000000000003697 <+13975>:	mulps  0x50(%rcx),%xmm15
   0x00000000000036ab <+13995>:	mulps  0x50(%rcx),%xmm12
   0x00000000000036b0 <+14000>:	addps  %xmm11,%xmm15
   0x00000000000036b4 <+14004>:	movaps 0x890(%r8),%xmm11
   0x00000000000036bc <+14012>:	mulps  0x50(%rcx),%xmm11
   0x00000000000036e3 <+14051>:	addps  %xmm10,%xmm12
   0x0000000000003705 <+14085>:	addps  %xmm13,%xmm11
   0x000000000000373c <+14140>:	movaps 0x8d0(%r8),%xmm15
   0x0000000000003744 <+14148>:	mulps  0x50(%rcx),%xmm15
   0x0000000000003749 <+14153>:	addps  %xmm13,%xmm15

1093
1094	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
1095	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
1096	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003672 <+13938>:	movaps 0x820(%r8),%xmm13
   0x000000000000369c <+13980>:	mulps  0x60(%rcx),%xmm13
   0x00000000000036c1 <+14017>:	addps  %xmm15,%xmm13
   0x00000000000036c5 <+14021>:	movaps 0x860(%r8),%xmm15
   0x00000000000036cd <+14029>:	mulps  0x60(%rcx),%xmm15
   0x00000000000036f4 <+14068>:	addps  %xmm12,%xmm15
   0x000000000000371a <+14106>:	movaps 0x8a0(%r8),%xmm15
   0x0000000000003722 <+14114>:	mulps  0x60(%rcx),%xmm15
   0x0000000000003727 <+14119>:	addps  %xmm11,%xmm15
   0x000000000000372b <+14123>:	movaps 0x8e0(%r8),%xmm11
   0x0000000000003733 <+14131>:	mulps  0x60(%rcx),%xmm11
   0x000000000000375a <+14170>:	addps  %xmm15,%xmm11

1097
1098	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
1099	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
1100	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000367a <+13946>:	movaps 0x830(%r8),%xmm14
   0x00000000000036a1 <+13985>:	mulps  0x70(%rcx),%xmm14
   0x00000000000036d2 <+14034>:	addps  %xmm13,%xmm14
   0x00000000000036e7 <+14055>:	movaps 0x8b0(%r8),%xmm10
   0x00000000000036ef <+14063>:	mulps  0x70(%rcx),%xmm10
   0x00000000000036f8 <+14072>:	movaps 0x870(%r8),%xmm12
   0x0000000000003700 <+14080>:	mulps  0x70(%rcx),%xmm12
   0x0000000000003716 <+14102>:	addps  %xmm15,%xmm12
   0x0000000000003738 <+14136>:	addps  %xmm15,%xmm10
   0x000000000000374d <+14157>:	movaps 0x8f0(%r8),%xmm13
   0x0000000000003755 <+14165>:	mulps  0x70(%rcx),%xmm13
   0x000000000000376e <+14190>:	addps  %xmm11,%xmm13

1101	      }
1102
1103	      /* A(3,2)*B(2,2) = C(3,2). */
1104	      for (i = 0; i < 4; i++)
1105	      {
1106	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1107	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
1108	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000375e <+14174>:	movaps 0x900(%r8),%xmm15
   0x0000000000003766 <+14182>:	mulps  0x140(%rcx),%xmm15
   0x0000000000003782 <+14210>:	addps  %xmm14,%xmm15
   0x000000000000379a <+14234>:	movaps 0x940(%r8),%xmm15
   0x00000000000037a2 <+14242>:	mulps  0x140(%rcx),%xmm15
   0x00000000000037be <+14270>:	addps  %xmm12,%xmm15
   0x00000000000037fe <+14334>:	movaps 0x980(%r8),%xmm12
   0x0000000000003806 <+14342>:	mulps  0x140(%rcx),%xmm12
   0x0000000000003822 <+14370>:	addps  %xmm10,%xmm12
   0x000000000000384e <+14414>:	movaps 0x9c0(%r8),%xmm10
   0x0000000000003856 <+14422>:	mulps  0x140(%rcx),%xmm10
   0x0000000000003872 <+14450>:	addps  %xmm13,%xmm10

1109
1110	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1111	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
1112	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003786 <+14214>:	movaps 0x910(%r8),%xmm14
   0x000000000000378e <+14222>:	mulps  0x150(%rcx),%xmm14
   0x0000000000003796 <+14230>:	addps  %xmm15,%xmm14
   0x00000000000037c2 <+14274>:	movaps 0x950(%r8),%xmm12
   0x00000000000037ca <+14282>:	mulps  0x150(%rcx),%xmm12
   0x00000000000037e6 <+14310>:	addps  %xmm15,%xmm12
   0x0000000000003826 <+14374>:	movaps 0x990(%r8),%xmm10
   0x000000000000382e <+14382>:	mulps  0x150(%rcx),%xmm10
   0x0000000000003836 <+14390>:	addps  %xmm12,%xmm10
   0x0000000000003876 <+14454>:	movaps 0x9d0(%r8),%xmm13
   0x000000000000387e <+14462>:	mulps  0x150(%rcx),%xmm13
   0x0000000000003886 <+14470>:	addps  %xmm10,%xmm13

1113
1114	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1115	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
1116	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003772 <+14194>:	movaps 0x920(%r8),%xmm11
   0x000000000000377a <+14202>:	mulps  0x160(%rcx),%xmm11
   0x00000000000037aa <+14250>:	addps  %xmm14,%xmm11
   0x00000000000037d6 <+14294>:	movaps 0x960(%r8),%xmm11
   0x00000000000037de <+14302>:	mulps  0x160(%rcx),%xmm11
   0x00000000000037fa <+14330>:	addps  %xmm12,%xmm11
   0x0000000000003812 <+14354>:	movaps 0x9a0(%r8),%xmm11
   0x000000000000381a <+14362>:	mulps  0x160(%rcx),%xmm11
   0x000000000000384a <+14410>:	addps  %xmm10,%xmm11
   0x0000000000003862 <+14434>:	movaps 0x9e0(%r8),%xmm11
   0x000000000000386a <+14442>:	mulps  0x160(%rcx),%xmm11
   0x000000000000389a <+14490>:	addps  %xmm13,%xmm11

1117
1118	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1119	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
1120	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000037ae <+14254>:	movaps 0x930(%r8),%xmm14
   0x00000000000037b6 <+14262>:	mulps  0x170(%rcx),%xmm14
   0x00000000000037d2 <+14290>:	addps  %xmm11,%xmm14
   0x00000000000037ea <+14314>:	movaps 0x970(%r8),%xmm15
   0x00000000000037f2 <+14322>:	mulps  0x170(%rcx),%xmm15
   0x000000000000380e <+14350>:	addps  %xmm11,%xmm15
   0x000000000000383a <+14394>:	movaps 0x9b0(%r8),%xmm12
   0x0000000000003842 <+14402>:	mulps  0x170(%rcx),%xmm12
   0x000000000000385e <+14430>:	addps  %xmm11,%xmm12
   0x000000000000388a <+14474>:	movaps 0x9f0(%r8),%xmm10
   0x0000000000003892 <+14482>:	mulps  0x170(%rcx),%xmm10
   0x00000000000038ae <+14510>:	addps  %xmm11,%xmm10

1121	      }
1122
1123	      /* A(3,3)*B(3,2) = C(3,2). */
1124	      for (i = 0; i < 4; i++)
1125	      {
1126	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1127	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
1128	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000389e <+14494>:	movaps 0xa00(%r8),%xmm13
   0x00000000000038a6 <+14502>:	mulps  0x240(%rcx),%xmm13
   0x00000000000038c2 <+14530>:	addps  %xmm14,%xmm13
   0x00000000000038ee <+14574>:	movaps 0xa40(%r8),%xmm14
   0x00000000000038f6 <+14582>:	mulps  0x240(%rcx),%xmm14
   0x0000000000003912 <+14610>:	addps  %xmm15,%xmm14
   0x000000000000393e <+14654>:	movaps 0xa80(%r8),%xmm15
   0x0000000000003946 <+14662>:	mulps  0x240(%rcx),%xmm15
   0x0000000000003962 <+14690>:	addps  %xmm12,%xmm15
   0x000000000000397a <+14714>:	movaps 0xac0(%r8),%xmm15
   0x0000000000003982 <+14722>:	mulps  0x240(%rcx),%xmm15
   0x000000000000399e <+14750>:	addps  %xmm10,%xmm15

1129
1130	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1131	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
1132	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000038c6 <+14534>:	movaps 0xa10(%r8),%xmm14
   0x00000000000038ce <+14542>:	mulps  0x250(%rcx),%xmm14
   0x00000000000038d6 <+14550>:	addps  %xmm13,%xmm14
   0x0000000000003916 <+14614>:	movaps 0xa50(%r8),%xmm15
   0x000000000000391e <+14622>:	mulps  0x250(%rcx),%xmm15
   0x0000000000003926 <+14630>:	addps  %xmm14,%xmm15
   0x0000000000003966 <+14694>:	movaps 0xa90(%r8),%xmm12
   0x000000000000396e <+14702>:	mulps  0x250(%rcx),%xmm12
   0x0000000000003976 <+14710>:	addps  %xmm15,%xmm12
   0x00000000000039a2 <+14754>:	movaps 0xad0(%r8),%xmm10
   0x00000000000039aa <+14762>:	mulps  0x250(%rcx),%xmm10
   0x00000000000039c6 <+14790>:	addps  %xmm15,%xmm10

1133
1134	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1135	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
1136	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000038da <+14554>:	movaps 0xa20(%r8),%xmm13
   0x00000000000038e2 <+14562>:	mulps  0x260(%rcx),%xmm13
   0x00000000000038ea <+14570>:	addps  %xmm14,%xmm13
   0x0000000000003902 <+14594>:	movaps 0xa60(%r8),%xmm13
   0x000000000000390a <+14602>:	mulps  0x260(%rcx),%xmm13
   0x000000000000393a <+14650>:	addps  %xmm15,%xmm13
   0x0000000000003952 <+14674>:	movaps 0xaa0(%r8),%xmm13
   0x000000000000395a <+14682>:	mulps  0x260(%rcx),%xmm13
   0x000000000000398a <+14730>:	addps  %xmm12,%xmm13
   0x00000000000039b6 <+14774>:	movaps 0xae0(%r8),%xmm13
   0x00000000000039be <+14782>:	mulps  0x260(%rcx),%xmm13
   0x00000000000039da <+14810>:	addps  %xmm10,%xmm13

1137
1138	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1139	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
1140	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000038b2 <+14514>:	movaps 0xa30(%r8),%xmm11
   0x00000000000038ba <+14522>:	mulps  0x270(%rcx),%xmm11
   0x00000000000038fe <+14590>:	addps  %xmm13,%xmm11
   0x000000000000392a <+14634>:	movaps 0xa70(%r8),%xmm14
   0x0000000000003932 <+14642>:	mulps  0x270(%rcx),%xmm14
   0x000000000000394e <+14670>:	addps  %xmm13,%xmm14
   0x000000000000398e <+14734>:	movaps 0xab0(%r8),%xmm12
   0x0000000000003996 <+14742>:	mulps  0x270(%rcx),%xmm12
   0x00000000000039b2 <+14770>:	addps  %xmm13,%xmm12
   0x00000000000039de <+14814>:	movaps 0xaf0(%r8),%xmm10
   0x00000000000039e6 <+14822>:	mulps  0x270(%rcx),%xmm10
   0x0000000000003a02 <+14850>:	addps  %xmm13,%xmm10

1141	      }
1142
1143	      /* A(3,4)*B(4,2) = C(3,2). */
1144	      for (i = 0; i < 4; i++)
1145	      {
1146	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1147	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
1148	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000039ca <+14794>:	movaps 0xb00(%r8),%xmm15
   0x00000000000039d2 <+14802>:	mulps  0x340(%rcx),%xmm15
   0x00000000000039ee <+14830>:	addps  %xmm11,%xmm15
   0x0000000000003a2e <+14894>:	movaps 0xb40(%r8),%xmm11
   0x0000000000003a36 <+14902>:	mulps  0x340(%rcx),%xmm11
   0x0000000000003a52 <+14930>:	addps  %xmm14,%xmm11
   0x0000000000003a7e <+14974>:	movaps 0xb80(%r8),%xmm14
   0x0000000000003a86 <+14982>:	mulps  0x340(%rcx),%xmm14
   0x0000000000003aa2 <+15010>:	addps  %xmm12,%xmm14
   0x0000000000003ace <+15054>:	movaps 0xbc0(%r8),%xmm12
   0x0000000000003ad6 <+15062>:	mulps  0x340(%rcx),%xmm12
   0x0000000000003af2 <+15090>:	addps  %xmm10,%xmm12

1149
1150	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1151	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
1152	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000039f2 <+14834>:	movaps 0xb10(%r8),%xmm11
   0x00000000000039fa <+14842>:	mulps  0x350(%rcx),%xmm11
   0x0000000000003a16 <+14870>:	addps  %xmm15,%xmm11
   0x0000000000003a56 <+14934>:	movaps 0xb50(%r8),%xmm14
   0x0000000000003a5e <+14942>:	mulps  0x350(%rcx),%xmm14
   0x0000000000003a66 <+14950>:	addps  %xmm11,%xmm14
   0x0000000000003aa6 <+15014>:	movaps 0xb90(%r8),%xmm12
   0x0000000000003aae <+15022>:	mulps  0x350(%rcx),%xmm12
   0x0000000000003ab6 <+15030>:	addps  %xmm14,%xmm12
   0x0000000000003af6 <+15094>:	movaps 0xbd0(%r8),%xmm10
   0x0000000000003afe <+15102>:	mulps  0x350(%rcx),%xmm10
   0x0000000000003b06 <+15110>:	addps  %xmm12,%xmm10

1153
1154	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1155	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
1156	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003a06 <+14854>:	movaps 0xb20(%r8),%xmm13
   0x0000000000003a0e <+14862>:	mulps  0x360(%rcx),%xmm13
   0x0000000000003a2a <+14890>:	addps  %xmm11,%xmm13
   0x0000000000003a6a <+14954>:	movaps 0xb60(%r8),%xmm11
   0x0000000000003a72 <+14962>:	mulps  0x360(%rcx),%xmm11
   0x0000000000003a7a <+14970>:	addps  %xmm14,%xmm11
   0x0000000000003aba <+15034>:	movaps 0xba0(%r8),%xmm14
   0x0000000000003ac2 <+15042>:	mulps  0x360(%rcx),%xmm14
   0x0000000000003aca <+15050>:	addps  %xmm12,%xmm14
   0x0000000000003b0a <+15114>:	movaps 0xbe0(%r8),%xmm12
   0x0000000000003b12 <+15122>:	mulps  0x360(%rcx),%xmm12
   0x0000000000003b1a <+15130>:	addps  %xmm10,%xmm12

1157
1158	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1159	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
1160	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003a1a <+14874>:	movaps 0xb30(%r8),%xmm15
   0x0000000000003a22 <+14882>:	mulps  0x370(%rcx),%xmm15
   0x0000000000003a3e <+14910>:	addps  %xmm13,%xmm15
   0x0000000000003a42 <+14914>:	movaps 0xb70(%r8),%xmm13
   0x0000000000003a4a <+14922>:	mulps  0x370(%rcx),%xmm13
   0x0000000000003a8e <+14990>:	addps  %xmm11,%xmm13
   0x0000000000003a92 <+14994>:	movaps 0xbb0(%r8),%xmm11
   0x0000000000003a9a <+15002>:	mulps  0x370(%rcx),%xmm11
   0x0000000000003ade <+15070>:	addps  %xmm14,%xmm11
   0x0000000000003ae2 <+15074>:	movaps 0xbf0(%r8),%xmm14
   0x0000000000003aea <+15082>:	mulps  0x370(%rcx),%xmm14
   0x0000000000003b60 <+15200>:	addps  %xmm12,%xmm14

1161	      }
1162
1163	      /* Store C(3,2) block. */
1164	      for (i = 0; i < 4; i++)
1165	      {
1166	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000003b1e <+15134>:	movaps 0x20(%rsp),%xmm10
   0x0000000000003b24 <+15140>:	mulps  %xmm10,%xmm15
   0x0000000000003b38 <+15160>:	mulps  %xmm10,%xmm13
   0x0000000000003b4c <+15180>:	mulps  %xmm10,%xmm11
   0x0000000000003b64 <+15204>:	mulps  %xmm10,%xmm14

1167	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_32]), C_row[i]);
   0x0000000000003b28 <+15144>:	addps  0x240(%r9),%xmm15
   0x0000000000003b3c <+15164>:	addps  0x250(%r9),%xmm13
   0x0000000000003b50 <+15184>:	addps  0x260(%r9),%xmm11
   0x0000000000003b68 <+15208>:	addps  0x270(%r9),%xmm14

1168	        _mm_store_ps(&C[i*4+C_OFFSET_32], C_row[i]);
   0x0000000000003b30 <+15152>:	movaps %xmm15,0x240(%r9)
   0x0000000000003b44 <+15172>:	movaps %xmm13,0x250(%r9)
   0x0000000000003b58 <+15192>:	movaps %xmm11,0x260(%r9)
   0x0000000000003b70 <+15216>:	movaps %xmm14,0x270(%r9)

1169	      }
1170	    }
1171
1172	    /* Reset C(3,3) matrix accumulators */
1173	    C_row[0] = _mm_setzero_ps();
1174	    C_row[1] = _mm_setzero_ps();
1175	    C_row[2] = _mm_setzero_ps();
1176	    C_row[3] = _mm_setzero_ps();
1177
1178	    if (norm_product[8][2] &&
   0x0000000000003b78 <+15224>:	comiss %xmm1,%xmm9
   0x0000000000003b7c <+15228>:	jb     0x40fe <stream_kernel+16638>

1179	        norm_product[9][6] &&
   0x0000000000003b82 <+15234>:	movss  0x90(%rsp),%xmm9
   0x0000000000003b8c <+15244>:	comiss %xmm1,%xmm9
   0x0000000000003b90 <+15248>:	jb     0x40fe <stream_kernel+16638>

1180	        norm_product[10][10] &&
   0x0000000000003b96 <+15254>:	movss  0x70(%rsp),%xmm9
   0x0000000000003b9d <+15261>:	comiss %xmm1,%xmm9
   0x0000000000003ba1 <+15265>:	jb     0x40fe <stream_kernel+16638>

1181	        norm_product[11][14])
   0x0000000000003ba7 <+15271>:	movss  0x60(%rsp),%xmm9
   0x0000000000003bae <+15278>:	comiss %xmm1,%xmm9
   0x0000000000003bb2 <+15282>:	jb     0x40fe <stream_kernel+16638>

1182	    {
1183	      /* A(3,1)*B(1,3) = C(3,3). */
1184	      for (i = 0; i < 4; i++)
1185	      {
1186	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
1187	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
1188	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003bb8 <+15288>:	movaps 0x800(%r8),%xmm9
   0x0000000000003bd8 <+15320>:	movaps 0x840(%r8),%xmm15
   0x0000000000003be8 <+15336>:	mulps  0x80(%rcx),%xmm9
   0x0000000000003c08 <+15368>:	mulps  0x80(%rcx),%xmm15
   0x0000000000003c38 <+15416>:	movaps 0x880(%r8),%xmm12
   0x0000000000003c48 <+15432>:	mulps  0x80(%rcx),%xmm12
   0x0000000000003c90 <+15504>:	movaps 0x8c0(%r8),%xmm12
   0x0000000000003c98 <+15512>:	mulps  0x80(%rcx),%xmm12

1189
1190	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
1191	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
1192	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003bc0 <+15296>:	movaps 0x810(%r8),%xmm12
   0x0000000000003be0 <+15328>:	movaps 0x850(%r8),%xmm11
   0x0000000000003bf0 <+15344>:	mulps  0x90(%rcx),%xmm12
   0x0000000000003c10 <+15376>:	mulps  0x90(%rcx),%xmm11
   0x0000000000003c20 <+15392>:	addps  %xmm9,%xmm12
   0x0000000000003c54 <+15444>:	movaps 0x890(%r8),%xmm13
   0x0000000000003c5c <+15452>:	mulps  0x90(%rcx),%xmm13
   0x0000000000003c64 <+15460>:	addps  %xmm15,%xmm11
   0x0000000000003c68 <+15464>:	movaps 0x8d0(%r8),%xmm15
   0x0000000000003c70 <+15472>:	mulps  0x90(%rcx),%xmm15
   0x0000000000003c8c <+15500>:	addps  %xmm12,%xmm13
   0x0000000000003cb4 <+15540>:	addps  %xmm12,%xmm15

1193
1194	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
1195	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
1196	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003bc8 <+15304>:	movaps 0x820(%r8),%xmm13
   0x0000000000003bf8 <+15352>:	mulps  0xa0(%rcx),%xmm13
   0x0000000000003c2c <+15404>:	movaps 0x860(%r8),%xmm9
   0x0000000000003c34 <+15412>:	addps  %xmm12,%xmm13
   0x0000000000003c40 <+15424>:	mulps  0xa0(%rcx),%xmm9
   0x0000000000003c78 <+15480>:	addps  %xmm11,%xmm9
   0x0000000000003ca4 <+15524>:	movaps 0x8a0(%r8),%xmm9
   0x0000000000003cac <+15532>:	mulps  0xa0(%rcx),%xmm9
   0x0000000000003cc8 <+15560>:	addps  %xmm13,%xmm9
   0x0000000000003ce0 <+15584>:	movaps 0x8e0(%r8),%xmm9
   0x0000000000003ce8 <+15592>:	mulps  0xa0(%rcx),%xmm9
   0x0000000000003d04 <+15620>:	addps  %xmm15,%xmm9

1197
1198	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
1199	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
1200	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003bd0 <+15312>:	movaps 0x830(%r8),%xmm10
   0x0000000000003c00 <+15360>:	mulps  0xb0(%rcx),%xmm10
   0x0000000000003c18 <+15384>:	movaps 0x8b0(%r8),%xmm14
   0x0000000000003c24 <+15396>:	mulps  0xb0(%rcx),%xmm14
   0x0000000000003c50 <+15440>:	addps  %xmm13,%xmm10
   0x0000000000003c7c <+15484>:	movaps 0x870(%r8),%xmm11
   0x0000000000003c84 <+15492>:	mulps  0xb0(%rcx),%xmm11
   0x0000000000003ca0 <+15520>:	addps  %xmm9,%xmm11
   0x0000000000003ccc <+15564>:	movaps 0x8f0(%r8),%xmm13
   0x0000000000003cd4 <+15572>:	mulps  0xb0(%rcx),%xmm13
   0x0000000000003cdc <+15580>:	addps  %xmm9,%xmm14
   0x0000000000003d18 <+15640>:	addps  %xmm9,%xmm13

1201	      }
1202
1203	      /* A(3,2)*B(2,3) = C(3,3). */
1204	      for (i = 0; i < 4; i++)
1205	      {
1206	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1207	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
1208	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003cb8 <+15544>:	movaps 0x900(%r8),%xmm12
   0x0000000000003cc0 <+15552>:	mulps  0x180(%rcx),%xmm12
   0x0000000000003cf0 <+15600>:	addps  %xmm10,%xmm12
   0x0000000000003d08 <+15624>:	movaps 0x940(%r8),%xmm15
   0x0000000000003d10 <+15632>:	mulps  0x180(%rcx),%xmm15
   0x0000000000003d68 <+15720>:	addps  %xmm11,%xmm15
   0x0000000000003d80 <+15744>:	movaps 0x980(%r8),%xmm15
   0x0000000000003d88 <+15752>:	mulps  0x180(%rcx),%xmm15
   0x0000000000003db8 <+15800>:	addps  %xmm14,%xmm15
   0x0000000000003dd0 <+15824>:	movaps 0x9c0(%r8),%xmm15
   0x0000000000003dd8 <+15832>:	mulps  0x180(%rcx),%xmm15
   0x0000000000003e08 <+15880>:	addps  %xmm13,%xmm15

1209
1210	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1211	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
1212	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003cf4 <+15604>:	movaps 0x910(%r8),%xmm10
   0x0000000000003cfc <+15612>:	mulps  0x190(%rcx),%xmm10
   0x0000000000003d2c <+15660>:	addps  %xmm12,%xmm10
   0x0000000000003d6c <+15724>:	movaps 0x950(%r8),%xmm11
   0x0000000000003d74 <+15732>:	mulps  0x190(%rcx),%xmm11
   0x0000000000003d7c <+15740>:	addps  %xmm15,%xmm11
   0x0000000000003dbc <+15804>:	movaps 0x990(%r8),%xmm14
   0x0000000000003dc4 <+15812>:	mulps  0x190(%rcx),%xmm14
   0x0000000000003dcc <+15820>:	addps  %xmm15,%xmm14
   0x0000000000003e0c <+15884>:	movaps 0x9d0(%r8),%xmm13
   0x0000000000003e14 <+15892>:	mulps  0x190(%rcx),%xmm13
   0x0000000000003e1c <+15900>:	addps  %xmm15,%xmm13

1213
1214	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1215	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
1216	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003d1c <+15644>:	movaps 0x920(%r8),%xmm9
   0x0000000000003d24 <+15652>:	mulps  0x1a0(%rcx),%xmm9
   0x0000000000003d40 <+15680>:	addps  %xmm10,%xmm9
   0x0000000000003d44 <+15684>:	movaps 0x960(%r8),%xmm10
   0x0000000000003d4c <+15692>:	mulps  0x1a0(%rcx),%xmm10
   0x0000000000003d90 <+15760>:	addps  %xmm11,%xmm10
   0x0000000000003d94 <+15764>:	movaps 0x9a0(%r8),%xmm11
   0x0000000000003d9c <+15772>:	mulps  0x1a0(%rcx),%xmm11
   0x0000000000003de0 <+15840>:	addps  %xmm14,%xmm11
   0x0000000000003de4 <+15844>:	movaps 0x9e0(%r8),%xmm14
   0x0000000000003dec <+15852>:	mulps  0x1a0(%rcx),%xmm14
   0x0000000000003e30 <+15920>:	addps  %xmm13,%xmm14

1217
1218	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1219	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
1220	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003d30 <+15664>:	movaps 0x930(%r8),%xmm12
   0x0000000000003d38 <+15672>:	mulps  0x1b0(%rcx),%xmm12
   0x0000000000003d54 <+15700>:	addps  %xmm9,%xmm12
   0x0000000000003d58 <+15704>:	movaps 0x970(%r8),%xmm9
   0x0000000000003d60 <+15712>:	mulps  0x1b0(%rcx),%xmm9
   0x0000000000003da4 <+15780>:	addps  %xmm10,%xmm9
   0x0000000000003da8 <+15784>:	movaps 0x9b0(%r8),%xmm10
   0x0000000000003db0 <+15792>:	mulps  0x1b0(%rcx),%xmm10
   0x0000000000003df4 <+15860>:	addps  %xmm11,%xmm10
   0x0000000000003df8 <+15864>:	movaps 0x9f0(%r8),%xmm11
   0x0000000000003e00 <+15872>:	mulps  0x1b0(%rcx),%xmm11
   0x0000000000003e44 <+15940>:	addps  %xmm14,%xmm11

1221	      }
1222
1223	      /* A(3,3)*B(3,3) = C(3,3). */
1224	      for (i = 0; i < 4; i++)
1225	      {
1226	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1227	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
1228	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003e20 <+15904>:	movaps 0xa00(%r8),%xmm15
   0x0000000000003e28 <+15912>:	mulps  0x280(%rcx),%xmm15
   0x0000000000003e48 <+15944>:	movaps 0xa40(%r8),%xmm14
   0x0000000000003e50 <+15952>:	mulps  0x280(%rcx),%xmm14
   0x0000000000003e58 <+15960>:	addps  %xmm12,%xmm15
   0x0000000000003e6c <+15980>:	addps  %xmm9,%xmm14
   0x0000000000003eac <+16044>:	movaps 0xa80(%r8),%xmm13
   0x0000000000003eb4 <+16052>:	mulps  0x280(%rcx),%xmm13
   0x0000000000003ee8 <+16104>:	movaps 0xac0(%r8),%xmm12
   0x0000000000003ef0 <+16112>:	mulps  0x280(%rcx),%xmm12
   0x0000000000003ef8 <+16120>:	addps  %xmm10,%xmm13
   0x0000000000003f0c <+16140>:	addps  %xmm11,%xmm12

1229
1230	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1231	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
1232	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003e5c <+15964>:	movaps 0xa10(%r8),%xmm12
   0x0000000000003e64 <+15972>:	mulps  0x290(%rcx),%xmm12
   0x0000000000003e70 <+15984>:	movaps 0xa50(%r8),%xmm9
   0x0000000000003e78 <+15992>:	mulps  0x290(%rcx),%xmm9
   0x0000000000003e80 <+16000>:	addps  %xmm15,%xmm12
   0x0000000000003ebc <+16060>:	addps  %xmm14,%xmm9
   0x0000000000003efc <+16124>:	movaps 0xa90(%r8),%xmm10
   0x0000000000003f04 <+16132>:	mulps  0x290(%rcx),%xmm10
   0x0000000000003f10 <+16144>:	movaps 0xad0(%r8),%xmm11
   0x0000000000003f18 <+16152>:	mulps  0x290(%rcx),%xmm11
   0x0000000000003f20 <+16160>:	addps  %xmm13,%xmm10
   0x0000000000003f5c <+16220>:	addps  %xmm12,%xmm11

1233
1234	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1235	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
1236	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003e34 <+15924>:	movaps 0xa20(%r8),%xmm13
   0x0000000000003e3c <+15932>:	mulps  0x2a0(%rcx),%xmm13
   0x0000000000003e94 <+16020>:	addps  %xmm12,%xmm13
   0x0000000000003e98 <+16024>:	movaps 0xa60(%r8),%xmm12
   0x0000000000003ea0 <+16032>:	mulps  0x2a0(%rcx),%xmm12
   0x0000000000003ed0 <+16080>:	addps  %xmm9,%xmm12
   0x0000000000003ed4 <+16084>:	movaps 0xaa0(%r8),%xmm9
   0x0000000000003edc <+16092>:	mulps  0x2a0(%rcx),%xmm9
   0x0000000000003f34 <+16180>:	addps  %xmm10,%xmm9
   0x0000000000003f38 <+16184>:	movaps 0xae0(%r8),%xmm10
   0x0000000000003f40 <+16192>:	mulps  0x2a0(%rcx),%xmm10
   0x0000000000003f70 <+16240>:	addps  %xmm11,%xmm10

1237
1238	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1239	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
1240	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003e84 <+16004>:	movaps 0xa30(%r8),%xmm15
   0x0000000000003e8c <+16012>:	mulps  0x2b0(%rcx),%xmm15
   0x0000000000003ea8 <+16040>:	addps  %xmm13,%xmm15
   0x0000000000003ec0 <+16064>:	movaps 0xa70(%r8),%xmm14
   0x0000000000003ec8 <+16072>:	mulps  0x2b0(%rcx),%xmm14
   0x0000000000003ee4 <+16100>:	addps  %xmm12,%xmm14
   0x0000000000003f24 <+16164>:	movaps 0xab0(%r8),%xmm13
   0x0000000000003f2c <+16172>:	mulps  0x2b0(%rcx),%xmm13
   0x0000000000003f48 <+16200>:	addps  %xmm9,%xmm13
   0x0000000000003f60 <+16224>:	movaps 0xaf0(%r8),%xmm12
   0x0000000000003f68 <+16232>:	mulps  0x2b0(%rcx),%xmm12
   0x0000000000003f84 <+16260>:	addps  %xmm10,%xmm12

1241	      }
1242
1243	      /* A(3,4)*B(4,3) = C(3,3). */
1244	      for (i = 0; i < 4; i++)
1245	      {
1246	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1247	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
1248	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003f4c <+16204>:	movaps 0xb00(%r8),%xmm9
   0x0000000000003f54 <+16212>:	mulps  0x380(%rcx),%xmm9
   0x0000000000003f98 <+16280>:	addps  %xmm15,%xmm9
   0x0000000000003fb0 <+16304>:	movaps 0xb40(%r8),%xmm9
   0x0000000000003fb8 <+16312>:	mulps  0x380(%rcx),%xmm9
   0x0000000000003fc4 <+16324>:	movaps 0xb80(%r8),%xmm15
   0x0000000000003fcc <+16332>:	mulps  0x380(%rcx),%xmm15
   0x0000000000003fe8 <+16360>:	addps  %xmm14,%xmm9
   0x0000000000003ffc <+16380>:	addps  %xmm13,%xmm15
   0x0000000000004028 <+16424>:	movaps 0xbc0(%r8),%xmm15
   0x0000000000004030 <+16432>:	mulps  0x380(%rcx),%xmm15
   0x0000000000004088 <+16520>:	addps  %xmm12,%xmm15

1249
1250	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1251	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
1252	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003f9c <+16284>:	movaps 0xb10(%r8),%xmm15
   0x0000000000003fa4 <+16292>:	mulps  0x390(%rcx),%xmm15
   0x0000000000003fac <+16300>:	addps  %xmm9,%xmm15
   0x0000000000003fec <+16364>:	movaps 0xb50(%r8),%xmm14
   0x0000000000003ff4 <+16372>:	mulps  0x390(%rcx),%xmm14
   0x0000000000004000 <+16384>:	movaps 0xb90(%r8),%xmm13
   0x0000000000004008 <+16392>:	mulps  0x390(%rcx),%xmm13
   0x0000000000004010 <+16400>:	addps  %xmm9,%xmm14
   0x0000000000004024 <+16420>:	addps  %xmm15,%xmm13
   0x000000000000408c <+16524>:	movaps 0xbd0(%r8),%xmm12
   0x0000000000004094 <+16532>:	mulps  0x390(%rcx),%xmm12
   0x000000000000409c <+16540>:	addps  %xmm15,%xmm12

1253
1254	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1255	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
1256	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003f88 <+16264>:	movaps 0xb20(%r8),%xmm10
   0x0000000000003f90 <+16272>:	mulps  0x3a0(%rcx),%xmm10
   0x0000000000003fc0 <+16320>:	addps  %xmm15,%xmm10
   0x0000000000004014 <+16404>:	movaps 0xb60(%r8),%xmm9
   0x000000000000401c <+16412>:	mulps  0x3a0(%rcx),%xmm9
   0x0000000000004038 <+16440>:	addps  %xmm14,%xmm9
   0x000000000000403c <+16444>:	movaps 0xba0(%r8),%xmm14
   0x0000000000004044 <+16452>:	mulps  0x3a0(%rcx),%xmm14
   0x0000000000004060 <+16480>:	addps  %xmm13,%xmm14
   0x0000000000004064 <+16484>:	movaps 0xbe0(%r8),%xmm13
   0x000000000000406c <+16492>:	mulps  0x3a0(%rcx),%xmm13
   0x00000000000040a0 <+16544>:	addps  %xmm12,%xmm13

1257
1258	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1259	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
1260	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003f74 <+16244>:	movaps 0xb30(%r8),%xmm11
   0x0000000000003f7c <+16252>:	mulps  0x3b0(%rcx),%xmm11
   0x0000000000003fd4 <+16340>:	addps  %xmm10,%xmm11
   0x0000000000003fd8 <+16344>:	movaps 0xb70(%r8),%xmm10
   0x0000000000003fe0 <+16352>:	mulps  0x3b0(%rcx),%xmm10
   0x000000000000404c <+16460>:	addps  %xmm9,%xmm10
   0x0000000000004050 <+16464>:	movaps 0xbb0(%r8),%xmm9
   0x0000000000004058 <+16472>:	mulps  0x3b0(%rcx),%xmm9
   0x0000000000004074 <+16500>:	addps  %xmm14,%xmm9
   0x0000000000004078 <+16504>:	movaps 0xbf0(%r8),%xmm14
   0x0000000000004080 <+16512>:	mulps  0x3b0(%rcx),%xmm14
   0x00000000000040d6 <+16598>:	addps  %xmm13,%xmm14

1261	      }
1262
1263	      /* Store C(3,3) block. */
1264	      for (i = 0; i < 4; i++)
1265	      {
1266	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000040a4 <+16548>:	movaps 0x20(%rsp),%xmm12
   0x00000000000040aa <+16554>:	mulps  %xmm12,%xmm11
   0x00000000000040b6 <+16566>:	mulps  %xmm12,%xmm10
   0x00000000000040c2 <+16578>:	mulps  %xmm12,%xmm9
   0x00000000000040da <+16602>:	mulps  %xmm12,%xmm14

1267	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_33]), C_row[i]);
   0x00000000000040ae <+16558>:	addps  0x280(%r9),%xmm11
   0x00000000000040ba <+16570>:	addps  0x290(%r9),%xmm10
   0x00000000000040c6 <+16582>:	addps  0x2a0(%r9),%xmm9
   0x00000000000040de <+16606>:	addps  0x2b0(%r9),%xmm14

1268	        _mm_store_ps(&C[i*4+C_OFFSET_33], C_row[i]);
   0x00000000000040ce <+16590>:	movaps %xmm11,0x280(%r9)
   0x00000000000040e6 <+16614>:	movaps %xmm10,0x290(%r9)
   0x00000000000040ee <+16622>:	movaps %xmm9,0x2a0(%r9)
   0x00000000000040f6 <+16630>:	movaps %xmm14,0x2b0(%r9)

1269	      }
1270	    }
1271
1272	    /* Reset C(3,4) matrix accumulators */
1273	    C_row[0] = _mm_setzero_ps();
1274	    C_row[1] = _mm_setzero_ps();
1275	    C_row[2] = _mm_setzero_ps();
1276	    C_row[3] = _mm_setzero_ps();
1277
1278	    if (norm_product[8][3] &&
   0x00000000000040fe <+16638>:	comiss %xmm1,%xmm8
   0x0000000000004102 <+16642>:	jb     0x4681 <stream_kernel+18049>

1279	        norm_product[9][7] &&
   0x0000000000004108 <+16648>:	movss  0x78(%rsp),%xmm8
   0x000000000000410f <+16655>:	comiss %xmm1,%xmm8
   0x0000000000004113 <+16659>:	jb     0x4681 <stream_kernel+18049>

1280	        norm_product[10][11] &&
   0x0000000000004119 <+16665>:	movss  0x68(%rsp),%xmm8
   0x0000000000004120 <+16672>:	comiss %xmm1,%xmm8
   0x0000000000004124 <+16676>:	jb     0x4681 <stream_kernel+18049>

1281	        norm_product[11][15])
   0x000000000000412a <+16682>:	movss  0x40(%rsp),%xmm8
   0x0000000000004131 <+16689>:	comiss %xmm1,%xmm8
   0x0000000000004135 <+16693>:	jb     0x4681 <stream_kernel+18049>

1282	    {
1283	      /* A(3,1)*B(1,4) = C(3,4). */
1284	      for (i = 0; i < 4; i++)
1285	      {
1286	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
1287	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
1288	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000413b <+16699>:	movaps 0x800(%r8),%xmm10
   0x000000000000415b <+16731>:	movaps 0x840(%r8),%xmm11
   0x000000000000417b <+16763>:	mulps  0xc0(%rcx),%xmm10
   0x000000000000419b <+16795>:	mulps  0xc0(%rcx),%xmm11
   0x00000000000041bf <+16831>:	movaps 0x880(%r8),%xmm10
   0x00000000000041c7 <+16839>:	mulps  0xc0(%rcx),%xmm10
   0x00000000000041fb <+16891>:	movaps 0x8c0(%r8),%xmm11
   0x0000000000004203 <+16899>:	mulps  0xc0(%rcx),%xmm11

1289
1290	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
1291	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
1292	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004143 <+16707>:	movaps 0x810(%r8),%xmm9
   0x0000000000004163 <+16739>:	movaps 0x850(%r8),%xmm12
   0x0000000000004183 <+16771>:	mulps  0xd0(%rcx),%xmm9
   0x00000000000041a3 <+16803>:	mulps  0xd0(%rcx),%xmm12
   0x00000000000041bb <+16827>:	addps  %xmm10,%xmm9
   0x00000000000041d3 <+16851>:	movaps 0x890(%r8),%xmm9
   0x00000000000041db <+16859>:	mulps  0xd0(%rcx),%xmm9
   0x00000000000041f7 <+16887>:	addps  %xmm11,%xmm12
   0x0000000000004233 <+16947>:	addps  %xmm10,%xmm9
   0x0000000000004237 <+16951>:	movaps 0x8d0(%r8),%xmm10
   0x000000000000423f <+16959>:	mulps  0xd0(%rcx),%xmm10
   0x000000000000426f <+17007>:	addps  %xmm11,%xmm10

1293
1294	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
1295	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
1296	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000414b <+16715>:	movaps 0x820(%r8),%xmm8
   0x000000000000416b <+16747>:	movaps 0x860(%r8),%xmm15
   0x000000000000418b <+16779>:	mulps  0xe0(%rcx),%xmm8
   0x00000000000041ab <+16811>:	mulps  0xe0(%rcx),%xmm15
   0x00000000000041cf <+16847>:	addps  %xmm9,%xmm8
   0x00000000000041e7 <+16871>:	movaps 0x8a0(%r8),%xmm8
   0x00000000000041ef <+16879>:	mulps  0xe0(%rcx),%xmm8
   0x000000000000420b <+16907>:	addps  %xmm12,%xmm15
   0x0000000000004247 <+16967>:	addps  %xmm9,%xmm8
   0x000000000000424b <+16971>:	movaps 0x8e0(%r8),%xmm9
   0x0000000000004253 <+16979>:	mulps  0xe0(%rcx),%xmm9
   0x0000000000004283 <+17027>:	addps  %xmm10,%xmm9

1297
1298	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
1299	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
1300	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004153 <+16723>:	movaps 0x830(%r8),%xmm14
   0x0000000000004173 <+16755>:	movaps 0x870(%r8),%xmm13
   0x0000000000004193 <+16787>:	mulps  0xf0(%rcx),%xmm14
   0x00000000000041b3 <+16819>:	mulps  0xf0(%rcx),%xmm13
   0x00000000000041e3 <+16867>:	addps  %xmm8,%xmm14
   0x000000000000420f <+16911>:	movaps 0x8b0(%r8),%xmm12
   0x0000000000004217 <+16919>:	mulps  0xf0(%rcx),%xmm12
   0x000000000000421f <+16927>:	addps  %xmm15,%xmm13
   0x000000000000425b <+16987>:	addps  %xmm8,%xmm12
   0x0000000000004273 <+17011>:	movaps 0x8f0(%r8),%xmm11
   0x000000000000427b <+17019>:	mulps  0xf0(%rcx),%xmm11
   0x0000000000004297 <+17047>:	addps  %xmm9,%xmm11

1301	      }
1302
1303	      /* A(3,2)*B(2,4) = C(3,4). */
1304	      for (i = 0; i < 4; i++)
1305	      {
1306	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1307	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
1308	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004223 <+16931>:	movaps 0x900(%r8),%xmm15
   0x000000000000422b <+16939>:	mulps  0x1c0(%rcx),%xmm15
   0x000000000000429b <+17051>:	movaps 0x940(%r8),%xmm9
   0x00000000000042a3 <+17059>:	mulps  0x1c0(%rcx),%xmm9
   0x00000000000042ab <+17067>:	addps  %xmm14,%xmm15
   0x00000000000042bf <+17087>:	addps  %xmm13,%xmm9
   0x00000000000042eb <+17131>:	movaps 0x980(%r8),%xmm14
   0x00000000000042f3 <+17139>:	mulps  0x1c0(%rcx),%xmm14
   0x000000000000434b <+17227>:	addps  %xmm12,%xmm14
   0x000000000000438b <+17291>:	movaps 0x9c0(%r8),%xmm12
   0x0000000000004393 <+17299>:	mulps  0x1c0(%rcx),%xmm12
   0x00000000000043af <+17327>:	addps  %xmm11,%xmm12

1309
1310	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1311	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
1312	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000042af <+17071>:	movaps 0x910(%r8),%xmm14
   0x00000000000042b7 <+17079>:	mulps  0x1d0(%rcx),%xmm14
   0x00000000000042c3 <+17091>:	movaps 0x950(%r8),%xmm13
   0x00000000000042cb <+17099>:	mulps  0x1d0(%rcx),%xmm13
   0x00000000000042d3 <+17107>:	addps  %xmm15,%xmm14
   0x000000000000430f <+17167>:	addps  %xmm9,%xmm13
   0x000000000000434f <+17231>:	movaps 0x990(%r8),%xmm12
   0x0000000000004357 <+17239>:	mulps  0x1d0(%rcx),%xmm12
   0x0000000000004373 <+17267>:	addps  %xmm14,%xmm12
   0x00000000000043b3 <+17331>:	movaps 0x9d0(%r8),%xmm11
   0x00000000000043bb <+17339>:	mulps  0x1d0(%rcx),%xmm11
   0x00000000000043d7 <+17367>:	addps  %xmm12,%xmm11

1313
1314	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1315	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
1316	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004287 <+17031>:	movaps 0x920(%r8),%xmm10
   0x000000000000428f <+17039>:	mulps  0x1e0(%rcx),%xmm10
   0x00000000000042d7 <+17111>:	movaps 0x9a0(%r8),%xmm15
   0x00000000000042df <+17119>:	mulps  0x1e0(%rcx),%xmm15
   0x00000000000042e7 <+17127>:	addps  %xmm14,%xmm10
   0x00000000000042ff <+17151>:	movaps 0x960(%r8),%xmm10
   0x0000000000004307 <+17159>:	mulps  0x1e0(%rcx),%xmm10
   0x0000000000004323 <+17187>:	addps  %xmm13,%xmm10
   0x0000000000004387 <+17287>:	addps  %xmm12,%xmm15
   0x00000000000043db <+17371>:	movaps 0x9e0(%r8),%xmm12
   0x00000000000043e3 <+17379>:	mulps  0x1e0(%rcx),%xmm12
   0x00000000000043eb <+17387>:	addps  %xmm11,%xmm12

1317
1318	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1319	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
1320	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000425f <+16991>:	movaps 0x930(%r8),%xmm8
   0x0000000000004267 <+16999>:	mulps  0x1f0(%rcx),%xmm8
   0x00000000000042fb <+17147>:	addps  %xmm10,%xmm8
   0x0000000000004313 <+17171>:	movaps 0x970(%r8),%xmm9
   0x000000000000431b <+17179>:	mulps  0x1f0(%rcx),%xmm9
   0x0000000000004337 <+17207>:	addps  %xmm10,%xmm9
   0x000000000000433b <+17211>:	movaps 0x9b0(%r8),%xmm10
   0x0000000000004343 <+17219>:	mulps  0x1f0(%rcx),%xmm10
   0x000000000000439b <+17307>:	addps  %xmm15,%xmm10
   0x000000000000439f <+17311>:	movaps 0x9f0(%r8),%xmm15
   0x00000000000043a7 <+17319>:	mulps  0x1f0(%rcx),%xmm15
   0x00000000000043ff <+17407>:	addps  %xmm12,%xmm15

1321	      }
1322
1323	      /* A(3,3)*B(3,4) = C(3,4). */
1324	      for (i = 0; i < 4; i++)
1325	      {
1326	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1327	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
1328	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004327 <+17191>:	movaps 0xa00(%r8),%xmm13
   0x000000000000432f <+17199>:	mulps  0x2c0(%rcx),%xmm13
   0x000000000000435f <+17247>:	addps  %xmm8,%xmm13
   0x0000000000004403 <+17411>:	movaps 0xa40(%r8),%xmm12
   0x000000000000440b <+17419>:	mulps  0x2c0(%rcx),%xmm12
   0x000000000000442b <+17451>:	movaps 0xa80(%r8),%xmm11
   0x0000000000004433 <+17459>:	mulps  0x2c0(%rcx),%xmm11
   0x000000000000443b <+17467>:	addps  %xmm9,%xmm12
   0x000000000000444f <+17487>:	addps  %xmm10,%xmm11
   0x000000000000448f <+17551>:	movaps 0xac0(%r8),%xmm8
   0x0000000000004497 <+17559>:	mulps  0x2c0(%rcx),%xmm8
   0x00000000000044db <+17627>:	addps  %xmm15,%xmm8

1329
1330	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1331	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
1332	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004363 <+17251>:	movaps 0xa10(%r8),%xmm8
   0x000000000000436b <+17259>:	mulps  0x2d0(%rcx),%xmm8
   0x00000000000043c3 <+17347>:	addps  %xmm13,%xmm8
   0x000000000000443f <+17471>:	movaps 0xa50(%r8),%xmm9
   0x0000000000004447 <+17479>:	mulps  0x2d0(%rcx),%xmm9
   0x0000000000004453 <+17491>:	movaps 0xa90(%r8),%xmm10
   0x000000000000445b <+17499>:	mulps  0x2d0(%rcx),%xmm10
   0x0000000000004463 <+17507>:	addps  %xmm12,%xmm9
   0x000000000000449f <+17567>:	addps  %xmm11,%xmm10
   0x00000000000044df <+17631>:	movaps 0xad0(%r8),%xmm15
   0x00000000000044e7 <+17639>:	mulps  0x2d0(%rcx),%xmm15
   0x0000000000004503 <+17667>:	addps  %xmm8,%xmm15

1333
1334	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1335	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
1336	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000043ef <+17391>:	movaps 0xa20(%r8),%xmm11
   0x00000000000043f7 <+17399>:	mulps  0x2e0(%rcx),%xmm11
   0x0000000000004413 <+17427>:	addps  %xmm8,%xmm11
   0x0000000000004417 <+17431>:	movaps 0xa60(%r8),%xmm8
   0x000000000000441f <+17439>:	mulps  0x2e0(%rcx),%xmm8
   0x0000000000004477 <+17527>:	addps  %xmm9,%xmm8
   0x000000000000447b <+17531>:	movaps 0xaa0(%r8),%xmm9
   0x0000000000004483 <+17539>:	mulps  0x2e0(%rcx),%xmm9
   0x00000000000044b3 <+17587>:	addps  %xmm10,%xmm9
   0x00000000000044b7 <+17591>:	movaps 0xae0(%r8),%xmm10
   0x00000000000044bf <+17599>:	mulps  0x2e0(%rcx),%xmm10
   0x0000000000004517 <+17687>:	addps  %xmm15,%xmm10

1337
1338	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1339	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
1340	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004377 <+17271>:	movaps 0xa30(%r8),%xmm14
   0x000000000000437f <+17279>:	mulps  0x2f0(%rcx),%xmm14
   0x00000000000043c7 <+17351>:	movaps 0xa70(%r8),%xmm13
   0x00000000000043cf <+17359>:	mulps  0x2f0(%rcx),%xmm13
   0x0000000000004427 <+17447>:	addps  %xmm11,%xmm14
   0x0000000000004467 <+17511>:	movaps 0xab0(%r8),%xmm12
   0x000000000000446f <+17519>:	mulps  0x2f0(%rcx),%xmm12
   0x000000000000448b <+17547>:	addps  %xmm8,%xmm13
   0x00000000000044a3 <+17571>:	movaps 0xaf0(%r8),%xmm11
   0x00000000000044ab <+17579>:	mulps  0x2f0(%rcx),%xmm11
   0x00000000000044c7 <+17607>:	addps  %xmm9,%xmm12
   0x000000000000452b <+17707>:	addps  %xmm10,%xmm11

1341	      }
1342
1343	      /* A(3,4)*B(4,4) = C(3,4). */
1344	      for (i = 0; i < 4; i++)
1345	      {
1346	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1347	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
1348	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000044cb <+17611>:	movaps 0xb00(%r8),%xmm9
   0x00000000000044d3 <+17619>:	mulps  0x3c0(%rcx),%xmm9
   0x00000000000044ef <+17647>:	addps  %xmm14,%xmm9
   0x0000000000004507 <+17671>:	movaps 0xb40(%r8),%xmm8
   0x000000000000450f <+17679>:	mulps  0x3c0(%rcx),%xmm8
   0x0000000000004553 <+17747>:	addps  %xmm13,%xmm8
   0x000000000000456b <+17771>:	movaps 0xb80(%r8),%xmm14
   0x0000000000004573 <+17779>:	mulps  0x3c0(%rcx),%xmm14
   0x00000000000045a3 <+17827>:	addps  %xmm12,%xmm14
   0x00000000000045e7 <+17895>:	movaps 0xbc0(%r8),%xmm12
   0x00000000000045ef <+17903>:	mulps  0x3c0(%rcx),%xmm12
   0x00000000000045fb <+17915>:	addps  %xmm11,%xmm12

1349
1350	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1351	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
1352	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000044f3 <+17651>:	movaps 0xb10(%r8),%xmm14
   0x00000000000044fb <+17659>:	mulps  0x3d0(%rcx),%xmm14
   0x000000000000453f <+17727>:	addps  %xmm9,%xmm14
   0x0000000000004557 <+17751>:	movaps 0xb50(%r8),%xmm13
   0x000000000000455f <+17759>:	mulps  0x3d0(%rcx),%xmm13
   0x000000000000458f <+17807>:	addps  %xmm8,%xmm13
   0x00000000000045a7 <+17831>:	movaps 0xb90(%r8),%xmm12
   0x00000000000045af <+17839>:	mulps  0x3d0(%rcx),%xmm12
   0x00000000000045df <+17887>:	addps  %xmm14,%xmm12
   0x00000000000045ff <+17919>:	movaps 0xbd0(%r8),%xmm11
   0x0000000000004607 <+17927>:	mulps  0x3d0(%rcx),%xmm11
   0x000000000000460f <+17935>:	addps  %xmm12,%xmm11

1353
1354	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1355	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
1356	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000451b <+17691>:	movaps 0xba0(%r8),%xmm15
   0x0000000000004523 <+17699>:	mulps  0x3e0(%rcx),%xmm15
   0x0000000000004543 <+17731>:	movaps 0xb20(%r8),%xmm9
   0x000000000000454b <+17739>:	mulps  0x3e0(%rcx),%xmm9
   0x0000000000004567 <+17767>:	addps  %xmm14,%xmm9
   0x0000000000004593 <+17811>:	movaps 0xb60(%r8),%xmm8
   0x000000000000459b <+17819>:	mulps  0x3e0(%rcx),%xmm8
   0x00000000000045b7 <+17847>:	addps  %xmm13,%xmm8
   0x00000000000045e3 <+17891>:	addps  %xmm12,%xmm15
   0x0000000000004613 <+17939>:	movaps 0xbe0(%r8),%xmm12
   0x000000000000461b <+17947>:	mulps  0x3e0(%rcx),%xmm12
   0x0000000000004623 <+17955>:	addps  %xmm11,%xmm12

1357
1358	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1359	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
1360	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000452f <+17711>:	movaps 0xb30(%r8),%xmm10
   0x0000000000004537 <+17719>:	mulps  0x3f0(%rcx),%xmm10
   0x000000000000457b <+17787>:	addps  %xmm9,%xmm10
   0x000000000000457f <+17791>:	movaps 0xb70(%r8),%xmm9
   0x0000000000004587 <+17799>:	mulps  0x3f0(%rcx),%xmm9
   0x00000000000045bb <+17851>:	movaps 0xbf0(%r8),%xmm13
   0x00000000000045c3 <+17859>:	mulps  0x3f0(%rcx),%xmm13
   0x00000000000045cb <+17867>:	addps  %xmm8,%xmm9
   0x00000000000045cf <+17871>:	movaps 0xbb0(%r8),%xmm8
   0x00000000000045d7 <+17879>:	mulps  0x3f0(%rcx),%xmm8
   0x00000000000045f7 <+17911>:	addps  %xmm15,%xmm8
   0x0000000000004659 <+18009>:	addps  %xmm12,%xmm13

1361	      }
1362
1363	      /* Store C(3,4) block. */
1364	      for (i = 0; i < 4; i++)
1365	      {
1366	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000004627 <+17959>:	movaps 0x20(%rsp),%xmm11
   0x000000000000462d <+17965>:	mulps  %xmm11,%xmm10
   0x0000000000004639 <+17977>:	mulps  %xmm11,%xmm9
   0x0000000000004645 <+17989>:	mulps  %xmm11,%xmm8
   0x000000000000465d <+18013>:	mulps  %xmm11,%xmm13

1367	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_34]), C_row[i]);
   0x0000000000004631 <+17969>:	addps  0x2c0(%r9),%xmm10
   0x000000000000463d <+17981>:	addps  0x2d0(%r9),%xmm9
   0x0000000000004649 <+17993>:	addps  0x2e0(%r9),%xmm8
   0x0000000000004661 <+18017>:	addps  0x2f0(%r9),%xmm13

1368	        _mm_store_ps(&C[i*4+C_OFFSET_34], C_row[i]);
   0x0000000000004651 <+18001>:	movaps %xmm10,0x2c0(%r9)
   0x0000000000004669 <+18025>:	movaps %xmm9,0x2d0(%r9)
   0x0000000000004671 <+18033>:	movaps %xmm8,0x2e0(%r9)
   0x0000000000004679 <+18041>:	movaps %xmm13,0x2f0(%r9)

1369	      }
1370	    }
1371
1372	    /* Reset C(4,1) matrix accumulators */
1373	    C_row[0] = _mm_setzero_ps();
1374	    C_row[1] = _mm_setzero_ps();
1375	    C_row[2] = _mm_setzero_ps();
1376	    C_row[3] = _mm_setzero_ps();
1377
1378	    if (norm_product[12][0] &&
   0x0000000000004681 <+18049>:	comiss %xmm1,%xmm7
   0x0000000000004684 <+18052>:	jb     0x4bc2 <stream_kernel+19394>

1379	        norm_product[13][4] &&
   0x000000000000468a <+18058>:	movss  0xb8(%rsp),%xmm7
   0x0000000000004693 <+18067>:	comiss %xmm1,%xmm7
   0x0000000000004696 <+18070>:	jb     0x4bc2 <stream_kernel+19394>

1380	        norm_product[14][8] &&
   0x000000000000469c <+18076>:	movss  0x58(%rsp),%xmm7
   0x00000000000046a2 <+18082>:	comiss %xmm1,%xmm7
   0x00000000000046a5 <+18085>:	jb     0x4bc2 <stream_kernel+19394>

1381	        norm_product[15][12])
   0x00000000000046ab <+18091>:	movss  0x38(%rsp),%xmm7
   0x00000000000046b1 <+18097>:	comiss %xmm1,%xmm7
   0x00000000000046b4 <+18100>:	jb     0x4bc2 <stream_kernel+19394>

1382	    {
1383	      /* A(4,1)*B(1,1) = C(4,1). */
1384	      for (i = 0; i < 4; i++)
1385	      {
1386	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1387	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
1388	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000046ba <+18106>:	movaps 0xc00(%r8),%xmm13
   0x00000000000046da <+18138>:	movaps 0xc40(%r8),%xmm8
   0x0000000000004702 <+18178>:	mulps  (%rcx),%xmm13
   0x0000000000004715 <+18197>:	mulps  (%rcx),%xmm8
   0x0000000000004730 <+18224>:	movaps 0xc80(%r8),%xmm13
   0x0000000000004738 <+18232>:	mulps  (%rcx),%xmm13
   0x0000000000004765 <+18277>:	movaps 0xcc0(%r8),%xmm8
   0x000000000000476d <+18285>:	mulps  (%rcx),%xmm8

1389
1390	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1391	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
1392	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000046c2 <+18114>:	movaps 0xc10(%r8),%xmm14
   0x00000000000046e2 <+18146>:	movaps 0xc50(%r8),%xmm7
   0x00000000000046fa <+18170>:	movaps 0xc90(%r8),%xmm15
   0x0000000000004706 <+18182>:	mulps  0x10(%rcx),%xmm14
   0x0000000000004719 <+18201>:	mulps  0x10(%rcx),%xmm7
   0x0000000000004727 <+18215>:	mulps  0x10(%rcx),%xmm15
   0x000000000000472c <+18220>:	addps  %xmm13,%xmm14
   0x0000000000004761 <+18273>:	addps  %xmm8,%xmm7
   0x0000000000004775 <+18293>:	movaps 0xcd0(%r8),%xmm7
   0x000000000000477d <+18301>:	mulps  0x10(%rcx),%xmm7
   0x0000000000004792 <+18322>:	addps  %xmm13,%xmm15
   0x00000000000047cb <+18379>:	addps  %xmm8,%xmm7

1393
1394	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1395	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
1396	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000046ca <+18122>:	movaps 0xc20(%r8),%xmm9
   0x00000000000046ea <+18154>:	movaps 0xc60(%r8),%xmm10
   0x000000000000470b <+18187>:	mulps  0x20(%rcx),%xmm9
   0x000000000000471d <+18205>:	mulps  0x20(%rcx),%xmm10
   0x000000000000473c <+18236>:	addps  %xmm14,%xmm9
   0x0000000000004754 <+18260>:	movaps 0xca0(%r8),%xmm9
   0x000000000000475c <+18268>:	mulps  0x20(%rcx),%xmm9
   0x0000000000004771 <+18289>:	addps  %xmm7,%xmm10
   0x0000000000004796 <+18326>:	movaps 0xce0(%r8),%xmm13
   0x000000000000479e <+18334>:	mulps  0x20(%rcx),%xmm13
   0x00000000000047a3 <+18339>:	addps  %xmm15,%xmm9
   0x00000000000047df <+18399>:	addps  %xmm7,%xmm13

1397
1398	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1399	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
1400	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000046d2 <+18130>:	movaps 0xc30(%r8),%xmm12
   0x00000000000046f2 <+18162>:	movaps 0xc70(%r8),%xmm11
   0x0000000000004710 <+18192>:	mulps  0x30(%rcx),%xmm12
   0x0000000000004722 <+18210>:	mulps  0x30(%rcx),%xmm11
   0x0000000000004750 <+18256>:	addps  %xmm9,%xmm12
   0x0000000000004781 <+18305>:	addps  %xmm10,%xmm11
   0x0000000000004785 <+18309>:	movaps 0xcb0(%r8),%xmm10
   0x000000000000478d <+18317>:	mulps  0x30(%rcx),%xmm10
   0x00000000000047b7 <+18359>:	addps  %xmm9,%xmm10
   0x00000000000047e3 <+18403>:	movaps 0xcf0(%r8),%xmm7
   0x00000000000047eb <+18411>:	mulps  0x30(%rcx),%xmm7
   0x0000000000004803 <+18435>:	addps  %xmm13,%xmm7

1401	      }
1402
1403	      /* A(4,2)*B(2,1) = C(4,1). */
1404	      for (i = 0; i < 4; i++)
1405	      {
1406	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1407	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
1408	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004740 <+18240>:	movaps 0xd00(%r8),%xmm14
   0x0000000000004748 <+18248>:	mulps  0x100(%rcx),%xmm14
   0x00000000000047ef <+18415>:	addps  %xmm12,%xmm14
   0x0000000000004807 <+18439>:	movaps 0xd40(%r8),%xmm13
   0x000000000000480f <+18447>:	mulps  0x100(%rcx),%xmm13
   0x0000000000004853 <+18515>:	addps  %xmm11,%xmm13
   0x000000000000486b <+18539>:	movaps 0xd80(%r8),%xmm13
   0x0000000000004873 <+18547>:	mulps  0x100(%rcx),%xmm13
   0x00000000000048a3 <+18595>:	addps  %xmm10,%xmm13
   0x00000000000048e3 <+18659>:	movaps 0xdc0(%r8),%xmm10
   0x00000000000048eb <+18667>:	mulps  0x100(%rcx),%xmm10
   0x0000000000004907 <+18695>:	addps  %xmm7,%xmm10

1409
1410	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1411	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
1412	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000047f3 <+18419>:	movaps 0xd10(%r8),%xmm12
   0x00000000000047fb <+18427>:	mulps  0x110(%rcx),%xmm12
   0x0000000000004817 <+18455>:	addps  %xmm14,%xmm12
   0x0000000000004857 <+18519>:	movaps 0xd50(%r8),%xmm11
   0x000000000000485f <+18527>:	mulps  0x110(%rcx),%xmm11
   0x0000000000004867 <+18535>:	addps  %xmm13,%xmm11
   0x00000000000048a7 <+18599>:	movaps 0xd90(%r8),%xmm10
   0x00000000000048af <+18607>:	mulps  0x110(%rcx),%xmm10
   0x00000000000048cb <+18635>:	addps  %xmm13,%xmm10
   0x000000000000490b <+18699>:	movaps 0xdd0(%r8),%xmm7
   0x0000000000004913 <+18707>:	mulps  0x110(%rcx),%xmm7
   0x000000000000492e <+18734>:	addps  %xmm10,%xmm7

1413
1414	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1415	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
1416	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000047bb <+18363>:	movaps 0xd20(%r8),%xmm9
   0x00000000000047c3 <+18371>:	mulps  0x120(%rcx),%xmm9
   0x000000000000481b <+18459>:	movaps 0xda0(%r8),%xmm14
   0x0000000000004823 <+18467>:	mulps  0x120(%rcx),%xmm14
   0x000000000000482b <+18475>:	addps  %xmm12,%xmm9
   0x000000000000482f <+18479>:	movaps 0xd60(%r8),%xmm12
   0x0000000000004837 <+18487>:	mulps  0x120(%rcx),%xmm12
   0x000000000000487b <+18555>:	addps  %xmm11,%xmm12
   0x00000000000048df <+18655>:	addps  %xmm10,%xmm14
   0x0000000000004932 <+18738>:	movaps 0xde0(%r8),%xmm10
   0x000000000000493a <+18746>:	mulps  0x120(%rcx),%xmm10
   0x0000000000004942 <+18754>:	addps  %xmm7,%xmm10

1417
1418	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1419	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
1420	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000047a7 <+18343>:	movaps 0xdb0(%r8),%xmm15
   0x00000000000047af <+18351>:	mulps  0x130(%rcx),%xmm15
   0x00000000000047cf <+18383>:	movaps 0xd30(%r8),%xmm8
   0x00000000000047d7 <+18391>:	mulps  0x130(%rcx),%xmm8
   0x000000000000483f <+18495>:	addps  %xmm9,%xmm8
   0x0000000000004843 <+18499>:	movaps 0xd70(%r8),%xmm9
   0x000000000000484b <+18507>:	mulps  0x130(%rcx),%xmm9
   0x000000000000488f <+18575>:	addps  %xmm12,%xmm9
   0x00000000000048f3 <+18675>:	addps  %xmm14,%xmm15
   0x00000000000048f7 <+18679>:	movaps 0xdf0(%r8),%xmm14
   0x00000000000048ff <+18687>:	mulps  0x130(%rcx),%xmm14
   0x0000000000004955 <+18773>:	addps  %xmm10,%xmm14

1421	      }
1422
1423	      /* A(4,3)*B(3,1) = C(4,1). */
1424	      for (i = 0; i < 4; i++)
1425	      {
1426	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1427	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
1428	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000487f <+18559>:	movaps 0xe00(%r8),%xmm11
   0x0000000000004887 <+18567>:	mulps  0x200(%rcx),%xmm11
   0x00000000000048b7 <+18615>:	addps  %xmm8,%xmm11
   0x0000000000004959 <+18777>:	movaps 0xe40(%r8),%xmm10
   0x0000000000004961 <+18785>:	mulps  0x200(%rcx),%xmm10
   0x0000000000004981 <+18817>:	movaps 0xe80(%r8),%xmm7
   0x0000000000004989 <+18825>:	mulps  0x200(%rcx),%xmm7
   0x0000000000004990 <+18832>:	addps  %xmm9,%xmm10
   0x00000000000049a4 <+18852>:	addps  %xmm15,%xmm7
   0x00000000000049e4 <+18916>:	movaps 0xec0(%r8),%xmm8
   0x00000000000049ec <+18924>:	mulps  0x200(%rcx),%xmm8
   0x0000000000004a25 <+18981>:	addps  %xmm14,%xmm8

1429
1430	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1431	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
1432	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000048bb <+18619>:	movaps 0xe10(%r8),%xmm8
   0x00000000000048c3 <+18627>:	mulps  0x210(%rcx),%xmm8
   0x000000000000491a <+18714>:	addps  %xmm11,%xmm8
   0x0000000000004994 <+18836>:	movaps 0xe50(%r8),%xmm9
   0x000000000000499c <+18844>:	mulps  0x210(%rcx),%xmm9
   0x00000000000049a8 <+18856>:	movaps 0xe90(%r8),%xmm15
   0x00000000000049b0 <+18864>:	mulps  0x210(%rcx),%xmm15
   0x00000000000049b8 <+18872>:	addps  %xmm10,%xmm9
   0x00000000000049f4 <+18932>:	addps  %xmm7,%xmm15
   0x00000000000049f8 <+18936>:	movaps 0xed0(%r8),%xmm7
   0x0000000000004a00 <+18944>:	mulps  0x210(%rcx),%xmm7
   0x0000000000004a39 <+19001>:	addps  %xmm8,%xmm7

1433
1434	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1435	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
1436	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004946 <+18758>:	movaps 0xe20(%r8),%xmm7
   0x000000000000494e <+18766>:	mulps  0x220(%rcx),%xmm7
   0x0000000000004969 <+18793>:	addps  %xmm8,%xmm7
   0x000000000000496d <+18797>:	movaps 0xe60(%r8),%xmm8
   0x0000000000004975 <+18805>:	mulps  0x220(%rcx),%xmm8
   0x00000000000049cc <+18892>:	addps  %xmm9,%xmm8
   0x00000000000049d0 <+18896>:	movaps 0xea0(%r8),%xmm9
   0x00000000000049d8 <+18904>:	mulps  0x220(%rcx),%xmm9
   0x0000000000004a07 <+18951>:	addps  %xmm15,%xmm9
   0x0000000000004a29 <+18985>:	movaps 0xee0(%r8),%xmm14
   0x0000000000004a31 <+18993>:	mulps  0x220(%rcx),%xmm14
   0x0000000000004a4d <+19021>:	addps  %xmm7,%xmm14

1437
1438	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1439	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
1440	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004893 <+18579>:	movaps 0xe70(%r8),%xmm12
   0x000000000000489b <+18587>:	mulps  0x230(%rcx),%xmm12
   0x00000000000048cf <+18639>:	movaps 0xe30(%r8),%xmm13
   0x00000000000048d7 <+18647>:	mulps  0x230(%rcx),%xmm13
   0x000000000000491e <+18718>:	movaps 0xeb0(%r8),%xmm11
   0x0000000000004926 <+18726>:	mulps  0x230(%rcx),%xmm11
   0x000000000000497d <+18813>:	addps  %xmm7,%xmm13
   0x00000000000049bc <+18876>:	movaps 0xef0(%r8),%xmm10
   0x00000000000049c4 <+18884>:	mulps  0x230(%rcx),%xmm10
   0x00000000000049e0 <+18912>:	addps  %xmm8,%xmm12
   0x0000000000004a11 <+18961>:	addps  %xmm9,%xmm11
   0x0000000000004a60 <+19040>:	addps  %xmm14,%xmm10

1441	      }
1442
1443	      /* A(4,4)*B(4,1) = C(4,1). */
1444	      for (i = 0; i < 4; i++)
1445	      {
1446	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1447	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
1448	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004a15 <+18965>:	movaps 0xf00(%r8),%xmm9
   0x0000000000004a1d <+18973>:	mulps  0x300(%rcx),%xmm9
   0x0000000000004a51 <+19025>:	movaps 0xf40(%r8),%xmm7
   0x0000000000004a59 <+19033>:	mulps  0x300(%rcx),%xmm7
   0x0000000000004a64 <+19044>:	addps  %xmm13,%xmm9
   0x0000000000004a78 <+19064>:	addps  %xmm12,%xmm7
   0x0000000000004aa4 <+19108>:	movaps 0xf80(%r8),%xmm13
   0x0000000000004aac <+19116>:	mulps  0x300(%rcx),%xmm13
   0x0000000000004aef <+19183>:	addps  %xmm11,%xmm13
   0x0000000000004b42 <+19266>:	movaps 0xfc0(%r8),%xmm13
   0x0000000000004b4a <+19274>:	mulps  0x300(%rcx),%xmm13
   0x0000000000004b8e <+19342>:	addps  %xmm10,%xmm13

1449
1450	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1451	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
1452	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004a68 <+19048>:	movaps 0xf10(%r8),%xmm13
   0x0000000000004a70 <+19056>:	mulps  0x310(%rcx),%xmm13
   0x0000000000004a7c <+19068>:	movaps 0xf50(%r8),%xmm12
   0x0000000000004a84 <+19076>:	mulps  0x310(%rcx),%xmm12
   0x0000000000004a8c <+19084>:	addps  %xmm9,%xmm13
   0x0000000000004adc <+19164>:	addps  %xmm7,%xmm12
   0x0000000000004af3 <+19187>:	movaps 0xf90(%r8),%xmm11
   0x0000000000004afb <+19195>:	mulps  0x310(%rcx),%xmm11
   0x0000000000004b3e <+19262>:	addps  %xmm13,%xmm11
   0x0000000000004b92 <+19346>:	movaps 0xfd0(%r8),%xmm10
   0x0000000000004b9a <+19354>:	mulps  0x310(%rcx),%xmm10
   0x0000000000004ba2 <+19362>:	addps  %xmm13,%xmm10

1453
1454	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1455	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
1456	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004a3d <+19005>:	movaps 0xf20(%r8),%xmm8
   0x0000000000004a45 <+19013>:	mulps  0x320(%rcx),%xmm8
   0x0000000000004aa0 <+19104>:	addps  %xmm13,%xmm8
   0x0000000000004ae0 <+19168>:	movaps 0xf60(%r8),%xmm7
   0x0000000000004ae8 <+19176>:	mulps  0x320(%rcx),%xmm7
   0x0000000000004b03 <+19203>:	addps  %xmm12,%xmm7
   0x0000000000004b07 <+19207>:	movaps 0xfa0(%r8),%xmm12
   0x0000000000004b0f <+19215>:	mulps  0x320(%rcx),%xmm12
   0x0000000000004b52 <+19282>:	addps  %xmm11,%xmm12
   0x0000000000004b56 <+19286>:	movaps 0xfe0(%r8),%xmm11
   0x0000000000004b5e <+19294>:	mulps  0x320(%rcx),%xmm11
   0x0000000000004ba6 <+19366>:	addps  %xmm10,%xmm11

1457
1458	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1459	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
1460	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004a90 <+19088>:	movaps 0xf30(%r8),%xmm9
   0x0000000000004a98 <+19096>:	mulps  0x330(%rcx),%xmm9
   0x0000000000004ab4 <+19124>:	addps  %xmm8,%xmm9
   0x0000000000004abc <+19132>:	movaps 0xf70(%r8),%xmm8
   0x0000000000004ac4 <+19140>:	mulps  0x330(%rcx),%xmm8
   0x0000000000004b17 <+19223>:	addps  %xmm7,%xmm8
   0x0000000000004b1f <+19231>:	movaps 0xfb0(%r8),%xmm7
   0x0000000000004b27 <+19239>:	mulps  0x330(%rcx),%xmm7
   0x0000000000004b66 <+19302>:	addps  %xmm12,%xmm7
   0x0000000000004b6e <+19310>:	movaps 0xff0(%r8),%xmm12
   0x0000000000004b76 <+19318>:	mulps  0x330(%rcx),%xmm12
   0x0000000000004baa <+19370>:	addps  %xmm11,%xmm12

1461	      }
1462
1463	      /* Store C(4,1) block. */
1464	      for (i = 0; i < 4; i++)
1465	      {
1466	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000004a0b <+18955>:	movaps 0x20(%rsp),%xmm15
   0x0000000000004ab8 <+19128>:	mulps  %xmm15,%xmm9
   0x0000000000004b1b <+19227>:	mulps  %xmm15,%xmm8
   0x0000000000004b6a <+19306>:	mulps  %xmm15,%xmm7
   0x0000000000004bae <+19374>:	mulps  %xmm15,%xmm12

1467	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_41]), C_row[i]);
   0x0000000000004acc <+19148>:	addps  0x300(%r9),%xmm9
   0x0000000000004b2e <+19246>:	addps  0x310(%r9),%xmm8
   0x0000000000004b7e <+19326>:	addps  0x320(%r9),%xmm7
   0x0000000000004bb2 <+19378>:	addps  0x330(%r9),%xmm12

1468	        _mm_store_ps(&C[i*4+C_OFFSET_41], C_row[i]);
   0x0000000000004ad4 <+19156>:	movaps %xmm9,0x300(%r9)
   0x0000000000004b36 <+19254>:	movaps %xmm8,0x310(%r9)
   0x0000000000004b86 <+19334>:	movaps %xmm7,0x320(%r9)
   0x0000000000004bba <+19386>:	movaps %xmm12,0x330(%r9)

1469	      }
1470	    }
1471
1472	    /* Reset C(4,2) matrix accumulators */
1473	    C_row[0] = _mm_setzero_ps();
1474	    C_row[1] = _mm_setzero_ps();
1475	    C_row[2] = _mm_setzero_ps();
1476	    C_row[3] = _mm_setzero_ps();
1477
1478	    if (norm_product[12][1] &&
   0x0000000000004bc2 <+19394>:	comiss %xmm1,%xmm6
   0x0000000000004bc5 <+19397>:	jb     0x5100 <stream_kernel+20736>

1479	        norm_product[13][5] &&
   0x0000000000004bcb <+19403>:	movss  0xa0(%rsp),%xmm6
   0x0000000000004bd4 <+19412>:	comiss %xmm1,%xmm6
   0x0000000000004bd7 <+19415>:	jb     0x5100 <stream_kernel+20736>

1480	        norm_product[14][9] &&
   0x0000000000004bdd <+19421>:	movss  0x50(%rsp),%xmm6
   0x0000000000004be3 <+19427>:	comiss %xmm1,%xmm6
   0x0000000000004be6 <+19430>:	jb     0x5100 <stream_kernel+20736>

1481	        norm_product[15][13])
   0x0000000000004bec <+19436>:	movss  0x18(%rsp),%xmm6
   0x0000000000004bf2 <+19442>:	comiss %xmm1,%xmm6
   0x0000000000004bf5 <+19445>:	jb     0x5100 <stream_kernel+20736>

1482	    {
1483	      /* A(4,1)*B(1,2) = C(4,2). */
1484	      for (i = 0; i < 4; i++)
1485	      {
1486	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1487	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
1488	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004bfb <+19451>:	movaps 0xc00(%r8),%xmm10
   0x0000000000004c1b <+19483>:	movaps 0xc40(%r8),%xmm13
   0x0000000000004c3b <+19515>:	movaps 0xc80(%r8),%xmm6
   0x0000000000004c43 <+19523>:	mulps  0x40(%rcx),%xmm10
   0x0000000000004c57 <+19543>:	mulps  0x40(%rcx),%xmm13
   0x0000000000004c6a <+19562>:	mulps  0x40(%rcx),%xmm6
   0x0000000000004ca1 <+19617>:	movaps 0xcc0(%r8),%xmm12
   0x0000000000004ca9 <+19625>:	mulps  0x40(%rcx),%xmm12

1489
1490	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1491	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
1492	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004c03 <+19459>:	movaps 0xc10(%r8),%xmm11
   0x0000000000004c23 <+19491>:	movaps 0xc50(%r8),%xmm14
   0x0000000000004c48 <+19528>:	mulps  0x50(%rcx),%xmm11
   0x0000000000004c5c <+19548>:	mulps  0x50(%rcx),%xmm14
   0x0000000000004c76 <+19574>:	addps  %xmm10,%xmm11
   0x0000000000004c7f <+19583>:	movaps 0xc90(%r8),%xmm10
   0x0000000000004c93 <+19603>:	mulps  0x50(%rcx),%xmm10
   0x0000000000004cae <+19630>:	addps  %xmm13,%xmm14
   0x0000000000004cb2 <+19634>:	movaps 0xcd0(%r8),%xmm13
   0x0000000000004cba <+19642>:	mulps  0x50(%rcx),%xmm13
   0x0000000000004ce3 <+19683>:	addps  %xmm6,%xmm10
   0x0000000000004d1b <+19739>:	addps  %xmm12,%xmm13

1493
1494	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1495	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
1496	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004c0b <+19467>:	movaps 0xc20(%r8),%xmm12
   0x0000000000004c2b <+19499>:	movaps 0xc60(%r8),%xmm7
   0x0000000000004c4d <+19533>:	mulps  0x60(%rcx),%xmm12
   0x0000000000004c61 <+19553>:	mulps  0x60(%rcx),%xmm7
   0x0000000000004c6e <+19566>:	movaps 0xce0(%r8),%xmm15
   0x0000000000004c7a <+19578>:	mulps  0x60(%rcx),%xmm15
   0x0000000000004c87 <+19591>:	addps  %xmm11,%xmm12
   0x0000000000004c8b <+19595>:	movaps 0xca0(%r8),%xmm11
   0x0000000000004c98 <+19608>:	mulps  0x60(%rcx),%xmm11
   0x0000000000004cbf <+19647>:	addps  %xmm14,%xmm7
   0x0000000000004cf3 <+19699>:	addps  %xmm10,%xmm11
   0x0000000000004d2f <+19759>:	addps  %xmm13,%xmm15

1497
1498	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1499	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
1500	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004c13 <+19475>:	movaps 0xc30(%r8),%xmm9
   0x0000000000004c33 <+19507>:	movaps 0xc70(%r8),%xmm8
   0x0000000000004c52 <+19538>:	mulps  0x70(%rcx),%xmm9
   0x0000000000004c65 <+19557>:	mulps  0x70(%rcx),%xmm8
   0x0000000000004c9d <+19613>:	addps  %xmm12,%xmm9
   0x0000000000004cd3 <+19667>:	addps  %xmm7,%xmm8
   0x0000000000004cd7 <+19671>:	movaps 0xcf0(%r8),%xmm7
   0x0000000000004cdf <+19679>:	mulps  0x70(%rcx),%xmm7
   0x0000000000004ce7 <+19687>:	movaps 0xcb0(%r8),%xmm6
   0x0000000000004cef <+19695>:	mulps  0x70(%rcx),%xmm6
   0x0000000000004d07 <+19719>:	addps  %xmm11,%xmm6
   0x0000000000004d43 <+19779>:	addps  %xmm15,%xmm7

1501	      }
1502
1503	      /* A(4,2)*B(2,2) = C(4,2). */
1504	      for (i = 0; i < 4; i++)
1505	      {
1506	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1507	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
1508	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004cf7 <+19703>:	movaps 0xd00(%r8),%xmm10
   0x0000000000004cff <+19711>:	mulps  0x140(%rcx),%xmm10
   0x0000000000004d0b <+19723>:	movaps 0xd40(%r8),%xmm11
   0x0000000000004d13 <+19731>:	mulps  0x140(%rcx),%xmm11
   0x0000000000004d1f <+19743>:	movaps 0xdc0(%r8),%xmm12
   0x0000000000004d27 <+19751>:	mulps  0x140(%rcx),%xmm12
   0x0000000000004d57 <+19799>:	addps  %xmm9,%xmm10
   0x0000000000004d6b <+19819>:	addps  %xmm8,%xmm11
   0x0000000000004dbf <+19903>:	movaps 0xd80(%r8),%xmm10
   0x0000000000004dc7 <+19911>:	mulps  0x140(%rcx),%xmm10
   0x0000000000004df7 <+19959>:	addps  %xmm6,%xmm10
   0x0000000000004e0a <+19978>:	addps  %xmm7,%xmm12

1509
1510	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1511	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
1512	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004d5b <+19803>:	movaps 0xd10(%r8),%xmm9
   0x0000000000004d63 <+19811>:	mulps  0x150(%rcx),%xmm9
   0x0000000000004d6f <+19823>:	movaps 0xd50(%r8),%xmm8
   0x0000000000004d77 <+19831>:	mulps  0x150(%rcx),%xmm8
   0x0000000000004d7f <+19839>:	addps  %xmm10,%xmm9
   0x0000000000004d93 <+19859>:	addps  %xmm11,%xmm8
   0x0000000000004dfb <+19963>:	movaps 0xd90(%r8),%xmm6
   0x0000000000004e03 <+19971>:	mulps  0x150(%rcx),%xmm6
   0x0000000000004e0e <+19982>:	movaps 0xdd0(%r8),%xmm7
   0x0000000000004e16 <+19990>:	mulps  0x150(%rcx),%xmm7
   0x0000000000004e1d <+19997>:	addps  %xmm10,%xmm6
   0x0000000000004e58 <+20056>:	addps  %xmm12,%xmm7

1513
1514	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1515	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
1516	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004d83 <+19843>:	movaps 0xd20(%r8),%xmm10
   0x0000000000004d8b <+19851>:	mulps  0x160(%rcx),%xmm10
   0x0000000000004d97 <+19863>:	movaps 0xda0(%r8),%xmm11
   0x0000000000004d9f <+19871>:	mulps  0x160(%rcx),%xmm11
   0x0000000000004da7 <+19879>:	addps  %xmm9,%xmm10
   0x0000000000004dab <+19883>:	movaps 0xd60(%r8),%xmm9
   0x0000000000004db3 <+19891>:	mulps  0x160(%rcx),%xmm9
   0x0000000000004dcf <+19919>:	addps  %xmm8,%xmm9
   0x0000000000004e31 <+20017>:	addps  %xmm6,%xmm11
   0x0000000000004e35 <+20021>:	movaps 0xde0(%r8),%xmm6
   0x0000000000004e3d <+20029>:	mulps  0x160(%rcx),%xmm6
   0x0000000000004e6c <+20076>:	addps  %xmm7,%xmm6

1517
1518	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1519	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
1520	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004cc3 <+19651>:	movaps 0xd70(%r8),%xmm14
   0x0000000000004ccb <+19659>:	mulps  0x170(%rcx),%xmm14
   0x0000000000004d33 <+19763>:	movaps 0xdb0(%r8),%xmm13
   0x0000000000004d3b <+19771>:	mulps  0x170(%rcx),%xmm13
   0x0000000000004d47 <+19783>:	movaps 0xd30(%r8),%xmm15
   0x0000000000004d4f <+19791>:	mulps  0x170(%rcx),%xmm15
   0x0000000000004dbb <+19899>:	addps  %xmm10,%xmm15
   0x0000000000004de3 <+19939>:	addps  %xmm9,%xmm14
   0x0000000000004e44 <+20036>:	addps  %xmm11,%xmm13
   0x0000000000004e5c <+20060>:	movaps 0xdf0(%r8),%xmm12
   0x0000000000004e64 <+20068>:	mulps  0x170(%rcx),%xmm12
   0x0000000000004e7e <+20094>:	addps  %xmm6,%xmm12

1521	      }
1522
1523	      /* A(4,3)*B(3,2) = C(4,2). */
1524	      for (i = 0; i < 4; i++)
1525	      {
1526	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1527	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
1528	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004dd3 <+19923>:	movaps 0xe00(%r8),%xmm8
   0x0000000000004ddb <+19931>:	mulps  0x240(%rcx),%xmm8
   0x0000000000004e82 <+20098>:	movaps 0xe40(%r8),%xmm6
   0x0000000000004e8a <+20106>:	mulps  0x240(%rcx),%xmm6
   0x0000000000004e91 <+20113>:	addps  %xmm15,%xmm8
   0x0000000000004ea5 <+20133>:	addps  %xmm14,%xmm6
   0x0000000000004ef8 <+20216>:	movaps 0xe80(%r8),%xmm6
   0x0000000000004f00 <+20224>:	mulps  0x240(%rcx),%xmm6
   0x0000000000004f2e <+20270>:	addps  %xmm13,%xmm6
   0x0000000000004f5a <+20314>:	movaps 0xec0(%r8),%xmm6
   0x0000000000004f62 <+20322>:	mulps  0x240(%rcx),%xmm6
   0x0000000000004f90 <+20368>:	addps  %xmm12,%xmm6

1529
1530	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1531	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
1532	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004e95 <+20117>:	movaps 0xe10(%r8),%xmm15
   0x0000000000004e9d <+20125>:	mulps  0x250(%rcx),%xmm15
   0x0000000000004ea9 <+20137>:	movaps 0xe50(%r8),%xmm14
   0x0000000000004eb1 <+20145>:	mulps  0x250(%rcx),%xmm14
   0x0000000000004eb9 <+20153>:	addps  %xmm8,%xmm15
   0x0000000000004ebd <+20157>:	movaps 0xed0(%r8),%xmm8
   0x0000000000004ec5 <+20165>:	mulps  0x250(%rcx),%xmm8
   0x0000000000004ef4 <+20212>:	addps  %xmm6,%xmm14
   0x0000000000004f32 <+20274>:	movaps 0xe90(%r8),%xmm13
   0x0000000000004f3a <+20282>:	mulps  0x250(%rcx),%xmm13
   0x0000000000004f56 <+20310>:	addps  %xmm6,%xmm13
   0x0000000000004fa4 <+20388>:	addps  %xmm6,%xmm8

1533
1534	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1535	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
1536	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004e6f <+20079>:	movaps 0xe20(%r8),%xmm7
   0x0000000000004e77 <+20087>:	mulps  0x260(%rcx),%xmm7
   0x0000000000004ecd <+20173>:	addps  %xmm15,%xmm7
   0x0000000000004ee5 <+20197>:	movaps 0xe60(%r8),%xmm7
   0x0000000000004eed <+20205>:	mulps  0x260(%rcx),%xmm7
   0x0000000000004f07 <+20231>:	addps  %xmm14,%xmm7
   0x0000000000004f1f <+20255>:	movaps 0xea0(%r8),%xmm7
   0x0000000000004f27 <+20263>:	mulps  0x260(%rcx),%xmm7
   0x0000000000004f69 <+20329>:	addps  %xmm13,%xmm7
   0x0000000000004f94 <+20372>:	movaps 0xee0(%r8),%xmm12
   0x0000000000004f9c <+20380>:	mulps  0x260(%rcx),%xmm12
   0x0000000000004fb7 <+20407>:	addps  %xmm8,%xmm12

1537
1538	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1539	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
1540	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004de7 <+19943>:	movaps 0xeb0(%r8),%xmm9
   0x0000000000004def <+19951>:	mulps  0x270(%rcx),%xmm9
   0x0000000000004e21 <+20001>:	movaps 0xe70(%r8),%xmm10
   0x0000000000004e29 <+20009>:	mulps  0x270(%rcx),%xmm10
   0x0000000000004e48 <+20040>:	movaps 0xe30(%r8),%xmm11
   0x0000000000004e50 <+20048>:	mulps  0x270(%rcx),%xmm11
   0x0000000000004ee1 <+20193>:	addps  %xmm7,%xmm11
   0x0000000000004f1b <+20251>:	addps  %xmm7,%xmm10
   0x0000000000004f7d <+20349>:	addps  %xmm7,%xmm9
   0x0000000000004fbb <+20411>:	movaps 0xef0(%r8),%xmm8
   0x0000000000004fc3 <+20419>:	mulps  0x270(%rcx),%xmm8
   0x0000000000004fdf <+20447>:	addps  %xmm12,%xmm8

1541	      }
1542
1543	      /* A(4,4)*B(4,2) = C(4,2). */
1544	      for (i = 0; i < 4; i++)
1545	      {
1546	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1547	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
1548	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004f0b <+20235>:	movaps 0xf00(%r8),%xmm14
   0x0000000000004f13 <+20243>:	mulps  0x340(%rcx),%xmm14
   0x0000000000004f42 <+20290>:	addps  %xmm11,%xmm14
   0x0000000000004fa8 <+20392>:	movaps 0xf40(%r8),%xmm6
   0x0000000000004fb0 <+20400>:	mulps  0x340(%rcx),%xmm6
   0x0000000000004fe3 <+20451>:	movaps 0xf80(%r8),%xmm12
   0x0000000000004feb <+20459>:	mulps  0x340(%rcx),%xmm12
   0x000000000000500b <+20491>:	addps  %xmm10,%xmm6
   0x000000000000501f <+20511>:	addps  %xmm9,%xmm12
   0x000000000000505e <+20574>:	movaps 0xfc0(%r8),%xmm11
   0x0000000000005066 <+20582>:	mulps  0x340(%rcx),%xmm11
   0x000000000000508a <+20618>:	addps  %xmm8,%xmm11

1549
1550	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1551	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
1552	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004f46 <+20294>:	movaps 0xf10(%r8),%xmm11
   0x0000000000004f4e <+20302>:	mulps  0x350(%rcx),%xmm11
   0x0000000000004fcb <+20427>:	addps  %xmm14,%xmm11
   0x000000000000500f <+20495>:	movaps 0xf50(%r8),%xmm10
   0x0000000000005017 <+20503>:	mulps  0x350(%rcx),%xmm10
   0x0000000000005023 <+20515>:	movaps 0xf90(%r8),%xmm9
   0x000000000000502b <+20523>:	mulps  0x350(%rcx),%xmm9
   0x0000000000005033 <+20531>:	addps  %xmm6,%xmm10
   0x000000000000506e <+20590>:	addps  %xmm12,%xmm9
   0x000000000000508e <+20622>:	movaps 0xfd0(%r8),%xmm8
   0x0000000000005096 <+20630>:	mulps  0x350(%rcx),%xmm8
   0x000000000000509e <+20638>:	addps  %xmm11,%xmm8

1553
1554	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1555	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
1556	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004ed1 <+20177>:	movaps 0xf20(%r8),%xmm15
   0x0000000000004ed9 <+20185>:	mulps  0x360(%rcx),%xmm15
   0x0000000000004ff3 <+20467>:	addps  %xmm11,%xmm15
   0x0000000000004ff7 <+20471>:	movaps 0xf60(%r8),%xmm11
   0x0000000000004fff <+20479>:	mulps  0x360(%rcx),%xmm11
   0x0000000000005046 <+20550>:	addps  %xmm10,%xmm11
   0x000000000000504a <+20554>:	movaps 0xfa0(%r8),%xmm10
   0x0000000000005052 <+20562>:	mulps  0x360(%rcx),%xmm10
   0x0000000000005072 <+20594>:	movaps 0xfe0(%r8),%xmm12
   0x000000000000507a <+20602>:	mulps  0x360(%rcx),%xmm12
   0x0000000000005082 <+20610>:	addps  %xmm9,%xmm10
   0x00000000000050a2 <+20642>:	addps  %xmm8,%xmm12

1557
1558	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1559	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
1560	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004f6d <+20333>:	movaps 0xfb0(%r8),%xmm13
   0x0000000000004f75 <+20341>:	mulps  0x370(%rcx),%xmm13
   0x0000000000004f81 <+20353>:	movaps 0xf30(%r8),%xmm7
   0x0000000000004f89 <+20361>:	mulps  0x370(%rcx),%xmm7
   0x0000000000004fcf <+20431>:	movaps 0xff0(%r8),%xmm14
   0x0000000000004fd7 <+20439>:	mulps  0x370(%rcx),%xmm14
   0x0000000000005007 <+20487>:	addps  %xmm15,%xmm7
   0x0000000000005037 <+20535>:	movaps 0xf70(%r8),%xmm6
   0x000000000000503f <+20543>:	mulps  0x370(%rcx),%xmm6
   0x000000000000505a <+20570>:	addps  %xmm11,%xmm6
   0x0000000000005086 <+20614>:	addps  %xmm10,%xmm13
   0x00000000000050a6 <+20646>:	addps  %xmm12,%xmm14

1561	      }
1562
1563	      /* Store C(4,2) block. */
1564	      for (i = 0; i < 4; i++)
1565	      {
1566	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000050aa <+20650>:	movaps 0x20(%rsp),%xmm12
   0x00000000000050b0 <+20656>:	mulps  %xmm12,%xmm7
   0x00000000000050bc <+20668>:	mulps  %xmm12,%xmm6
   0x00000000000050c8 <+20680>:	mulps  %xmm12,%xmm13
   0x00000000000050d4 <+20692>:	mulps  %xmm12,%xmm14

1567	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_42]), C_row[i]);
   0x00000000000050b4 <+20660>:	addps  0x340(%r9),%xmm7
   0x00000000000050c0 <+20672>:	addps  0x350(%r9),%xmm6
   0x00000000000050cc <+20684>:	addps  0x360(%r9),%xmm13
   0x00000000000050d8 <+20696>:	addps  0x370(%r9),%xmm14

1568	        _mm_store_ps(&C[i*4+C_OFFSET_42], C_row[i]);
   0x00000000000050e0 <+20704>:	movaps %xmm7,0x340(%r9)
   0x00000000000050e8 <+20712>:	movaps %xmm6,0x350(%r9)
   0x00000000000050f0 <+20720>:	movaps %xmm13,0x360(%r9)
   0x00000000000050f8 <+20728>:	movaps %xmm14,0x370(%r9)

1569	      }
1570	    }
1571
1572	    /* Reset C(4,3) matrix accumulators */
1573	    C_row[0] = _mm_setzero_ps();
1574	    C_row[1] = _mm_setzero_ps();
1575	    C_row[2] = _mm_setzero_ps();
1576	    C_row[3] = _mm_setzero_ps();
1577
1578	    if (norm_product[12][2] &&
   0x0000000000005100 <+20736>:	comiss %xmm1,%xmm5
   0x0000000000005103 <+20739>:	jb     0x55f5 <stream_kernel+22005>

1579	        norm_product[13][6] &&
   0x0000000000005109 <+20745>:	movss  0x80(%rsp),%xmm5
   0x0000000000005112 <+20754>:	comiss %xmm1,%xmm5
   0x0000000000005115 <+20757>:	jb     0x55f5 <stream_kernel+22005>

1580	        norm_product[14][10] &&
   0x000000000000511b <+20763>:	movss  0x48(%rsp),%xmm5
   0x0000000000005121 <+20769>:	comiss %xmm1,%xmm5
   0x0000000000005124 <+20772>:	jb     0x55f5 <stream_kernel+22005>

1581	        norm_product[15][14])
   0x000000000000512a <+20778>:	movss  0x30(%rsp),%xmm5
   0x0000000000005130 <+20784>:	comiss %xmm1,%xmm5
   0x0000000000005133 <+20787>:	jb     0x55f5 <stream_kernel+22005>

1582	    {
1583	      /* A(4,1)*B(1,3) = C(4,3). */
1584	      for (i = 0; i < 4; i++)
1585	      {
1586	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1587	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
1588	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005139 <+20793>:	movaps 0x80(%rcx),%xmm8
   0x0000000000005159 <+20825>:	movaps 0xc00(%r8),%xmm7
   0x0000000000005171 <+20849>:	movaps 0xc40(%r8),%xmm9
   0x0000000000005191 <+20881>:	mulps  %xmm8,%xmm7
   0x00000000000051bd <+20925>:	movaps 0xc80(%r8),%xmm14
   0x00000000000051c5 <+20933>:	mulps  %xmm8,%xmm9
   0x00000000000051f9 <+20985>:	mulps  %xmm8,%xmm14
   0x0000000000005205 <+20997>:	movaps 0xcc0(%r8),%xmm14
   0x000000000000522c <+21036>:	mulps  %xmm8,%xmm14

1589
1590	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1591	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
1592	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005141 <+20801>:	movaps 0x90(%rcx),%xmm13
   0x0000000000005161 <+20833>:	movaps 0xc10(%r8),%xmm12
   0x0000000000005179 <+20857>:	movaps 0xc50(%r8),%xmm5
   0x0000000000005195 <+20885>:	mulps  %xmm13,%xmm12
   0x0000000000005199 <+20889>:	addps  %xmm7,%xmm12
   0x00000000000051c9 <+20937>:	mulps  %xmm13,%xmm5
   0x00000000000051cd <+20941>:	addps  %xmm9,%xmm5
   0x00000000000051d1 <+20945>:	movaps 0xc90(%r8),%xmm9
   0x00000000000051fd <+20989>:	mulps  %xmm13,%xmm9
   0x0000000000005201 <+20993>:	addps  %xmm14,%xmm9
   0x0000000000005230 <+21040>:	movaps 0xcd0(%r8),%xmm8
   0x0000000000005238 <+21048>:	mulps  %xmm13,%xmm8
   0x0000000000005244 <+21060>:	addps  %xmm14,%xmm8

1593
1594	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1595	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
1596	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005149 <+20809>:	movaps 0xa0(%rcx),%xmm11
   0x0000000000005169 <+20841>:	movaps 0xc20(%r8),%xmm14
   0x0000000000005181 <+20865>:	movaps 0xc60(%r8),%xmm15
   0x00000000000051a5 <+20901>:	mulps  %xmm11,%xmm14
   0x00000000000051a9 <+20905>:	addps  %xmm12,%xmm14
   0x00000000000051d9 <+20953>:	mulps  %xmm11,%xmm15
   0x00000000000051dd <+20957>:	addps  %xmm5,%xmm15
   0x00000000000051e1 <+20961>:	movaps 0xca0(%r8),%xmm5
   0x000000000000520d <+21005>:	mulps  %xmm11,%xmm5
   0x0000000000005211 <+21009>:	addps  %xmm9,%xmm5
   0x000000000000523c <+21052>:	movaps 0xce0(%r8),%xmm13
   0x0000000000005248 <+21064>:	mulps  %xmm11,%xmm13
   0x0000000000005254 <+21076>:	addps  %xmm8,%xmm13

1597
1598	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1599	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
1600	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005151 <+20817>:	movaps 0xb0(%rcx),%xmm10
   0x0000000000005189 <+20873>:	movaps 0xc70(%r8),%xmm6
   0x000000000000519d <+20893>:	movaps 0xcb0(%r8),%xmm7
   0x00000000000051ad <+20909>:	movaps 0xc30(%r8),%xmm12
   0x00000000000051b5 <+20917>:	mulps  %xmm10,%xmm12
   0x00000000000051b9 <+20921>:	addps  %xmm14,%xmm12
   0x00000000000051e9 <+20969>:	mulps  %xmm10,%xmm6
   0x00000000000051ed <+20973>:	addps  %xmm15,%xmm6
   0x000000000000521d <+21021>:	mulps  %xmm10,%xmm7
   0x0000000000005221 <+21025>:	addps  %xmm5,%xmm7
   0x0000000000005224 <+21028>:	movaps 0xcf0(%r8),%xmm5
   0x0000000000005260 <+21088>:	mulps  %xmm10,%xmm5
   0x0000000000005270 <+21104>:	addps  %xmm13,%xmm5

1601	      }
1602
1603	      /* A(4,2)*B(2,3) = C(4,3). */
1604	      for (i = 0; i < 4; i++)
1605	      {
1606	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1607	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
1608	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005215 <+21013>:	movaps 0x180(%rcx),%xmm9
   0x000000000000524c <+21068>:	movaps 0xd00(%r8),%xmm11
   0x000000000000526c <+21100>:	mulps  %xmm9,%xmm11
   0x0000000000005280 <+21120>:	addps  %xmm12,%xmm11
   0x00000000000052a0 <+21152>:	movaps 0xd40(%r8),%xmm11
   0x00000000000052b4 <+21172>:	mulps  %xmm9,%xmm11
   0x00000000000052bc <+21180>:	movaps 0xd80(%r8),%xmm13
   0x00000000000052c8 <+21192>:	addps  %xmm6,%xmm11
   0x00000000000052ec <+21228>:	mulps  %xmm9,%xmm13
   0x00000000000052f0 <+21232>:	movaps 0xdc0(%r8),%xmm6
   0x0000000000005308 <+21256>:	addps  %xmm7,%xmm13
   0x000000000000532c <+21292>:	mulps  %xmm9,%xmm6
   0x0000000000005348 <+21320>:	addps  %xmm5,%xmm6

1609
1610	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1611	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
1612	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005264 <+21092>:	movaps 0xd10(%r8),%xmm10
   0x0000000000005284 <+21124>:	movaps 0x190(%rcx),%xmm14
   0x0000000000005294 <+21140>:	mulps  %xmm14,%xmm10
   0x0000000000005298 <+21144>:	addps  %xmm11,%xmm10
   0x00000000000052cc <+21196>:	movaps 0xd50(%r8),%xmm6
   0x00000000000052d4 <+21204>:	mulps  %xmm14,%xmm6
   0x00000000000052d8 <+21208>:	addps  %xmm11,%xmm6
   0x000000000000530c <+21260>:	movaps 0xd90(%r8),%xmm7
   0x0000000000005314 <+21268>:	mulps  %xmm14,%xmm7
   0x0000000000005318 <+21272>:	addps  %xmm13,%xmm7
   0x0000000000005330 <+21296>:	movaps 0xdd0(%r8),%xmm9
   0x000000000000533c <+21308>:	mulps  %xmm14,%xmm9
   0x0000000000005357 <+21335>:	addps  %xmm6,%xmm9

1613
1614	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1615	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
1616	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000051f1 <+20977>:	movaps 0x1a0(%rcx),%xmm15
   0x0000000000005274 <+21108>:	movaps 0xd20(%r8),%xmm13
   0x000000000000527c <+21116>:	mulps  %xmm15,%xmm13
   0x00000000000052a8 <+21160>:	addps  %xmm10,%xmm13
   0x00000000000052dc <+21212>:	movaps 0xd60(%r8),%xmm11
   0x00000000000052e4 <+21220>:	mulps  %xmm15,%xmm11
   0x00000000000052e8 <+21224>:	addps  %xmm6,%xmm11
   0x000000000000531c <+21276>:	movaps 0xda0(%r8),%xmm13
   0x0000000000005324 <+21284>:	mulps  %xmm15,%xmm13
   0x0000000000005328 <+21288>:	addps  %xmm7,%xmm13
   0x0000000000005340 <+21312>:	movaps 0xde0(%r8),%xmm14
   0x000000000000534b <+21323>:	mulps  %xmm15,%xmm14
   0x000000000000537c <+21372>:	addps  %xmm9,%xmm14

1617
1618	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1619	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
1620	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005258 <+21080>:	movaps 0xd30(%r8),%xmm8
   0x000000000000528c <+21132>:	movaps 0x1b0(%rcx),%xmm12
   0x000000000000529c <+21148>:	mulps  %xmm12,%xmm8
   0x00000000000052ac <+21164>:	movaps 0xd70(%r8),%xmm10
   0x00000000000052b8 <+21176>:	addps  %xmm13,%xmm8
   0x00000000000052c4 <+21188>:	mulps  %xmm12,%xmm10
   0x00000000000052f8 <+21240>:	addps  %xmm11,%xmm10
   0x00000000000052fc <+21244>:	movaps 0xdb0(%r8),%xmm11
   0x0000000000005304 <+21252>:	mulps  %xmm12,%xmm11
   0x0000000000005338 <+21304>:	addps  %xmm13,%xmm11
   0x0000000000005380 <+21376>:	movaps 0xdf0(%r8),%xmm9
   0x0000000000005388 <+21384>:	mulps  %xmm12,%xmm9
   0x0000000000005398 <+21400>:	addps  %xmm14,%xmm9

1621	      }
1622
1623	      /* A(4,3)*B(3,3) = C(4,3). */
1624	      for (i = 0; i < 4; i++)
1625	      {
1626	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1627	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
1628	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000535b <+21339>:	movaps 0x280(%rcx),%xmm5
   0x000000000000538c <+21388>:	movaps 0xe00(%r8),%xmm12
   0x0000000000005394 <+21396>:	mulps  %xmm5,%xmm12
   0x00000000000053a8 <+21416>:	addps  %xmm8,%xmm12
   0x00000000000053ac <+21420>:	movaps 0xe40(%r8),%xmm8
   0x00000000000053b4 <+21428>:	mulps  %xmm5,%xmm8
   0x00000000000053cc <+21452>:	movaps 0xe80(%r8),%xmm15
   0x00000000000053d4 <+21460>:	mulps  %xmm5,%xmm15
   0x00000000000053e8 <+21480>:	addps  %xmm10,%xmm8
   0x00000000000053fc <+21500>:	movaps 0xec0(%r8),%xmm8
   0x0000000000005404 <+21508>:	mulps  %xmm5,%xmm8
   0x0000000000005423 <+21539>:	addps  %xmm11,%xmm15
   0x000000000000548f <+21647>:	addps  %xmm9,%xmm8

1629
1630	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1631	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
1632	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000534f <+21327>:	movaps 0xe10(%r8),%xmm15
   0x0000000000005362 <+21346>:	movaps 0x290(%rcx),%xmm6
   0x0000000000005378 <+21368>:	mulps  %xmm6,%xmm15
   0x00000000000053b8 <+21432>:	addps  %xmm12,%xmm15
   0x00000000000053dc <+21468>:	movaps 0xe50(%r8),%xmm14
   0x00000000000053e4 <+21476>:	mulps  %xmm6,%xmm14
   0x00000000000053f8 <+21496>:	addps  %xmm8,%xmm14
   0x0000000000005408 <+21512>:	movaps 0xed0(%r8),%xmm5
   0x0000000000005410 <+21520>:	mulps  %xmm6,%xmm5
   0x0000000000005427 <+21543>:	movaps 0xe90(%r8),%xmm11
   0x000000000000542f <+21551>:	mulps  %xmm6,%xmm11
   0x0000000000005443 <+21571>:	addps  %xmm15,%xmm11
   0x000000000000549b <+21659>:	addps  %xmm8,%xmm5

1633
1634	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1635	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
1636	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005369 <+21353>:	movaps 0x2a0(%rcx),%xmm13
   0x000000000000539c <+21404>:	movaps 0xe20(%r8),%xmm14
   0x00000000000053a4 <+21412>:	mulps  %xmm13,%xmm14
   0x00000000000053c8 <+21448>:	addps  %xmm15,%xmm14
   0x00000000000053ec <+21484>:	movaps 0xe60(%r8),%xmm10
   0x00000000000053f4 <+21492>:	mulps  %xmm13,%xmm10
   0x0000000000005413 <+21523>:	addps  %xmm14,%xmm10
   0x0000000000005437 <+21559>:	movaps 0xea0(%r8),%xmm10
   0x000000000000543f <+21567>:	mulps  %xmm13,%xmm10
   0x0000000000005453 <+21587>:	addps  %xmm11,%xmm10
   0x0000000000005457 <+21591>:	movaps 0xee0(%r8),%xmm11
   0x0000000000005467 <+21607>:	mulps  %xmm13,%xmm11
   0x00000000000054ab <+21675>:	addps  %xmm5,%xmm11

1637
1638	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1639	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
1640	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005371 <+21361>:	movaps 0x2b0(%rcx),%xmm7
   0x00000000000053bc <+21436>:	movaps 0xe30(%r8),%xmm12
   0x00000000000053c4 <+21444>:	mulps  %xmm7,%xmm12
   0x00000000000053d8 <+21464>:	addps  %xmm14,%xmm12
   0x0000000000005417 <+21527>:	movaps 0xe70(%r8),%xmm14
   0x000000000000541f <+21535>:	mulps  %xmm7,%xmm14
   0x0000000000005433 <+21555>:	addps  %xmm10,%xmm14
   0x0000000000005447 <+21575>:	movaps 0xeb0(%r8),%xmm15
   0x000000000000544f <+21583>:	mulps  %xmm7,%xmm15
   0x000000000000546b <+21611>:	movaps 0xef0(%r8),%xmm13
   0x0000000000005473 <+21619>:	mulps  %xmm7,%xmm13
   0x000000000000547f <+21631>:	addps  %xmm10,%xmm15
   0x00000000000054ba <+21690>:	addps  %xmm11,%xmm13

1641	      }
1642
1643	      /* A(4,4)*B(4,3) = C(4,3). */
1644	      for (i = 0; i < 4; i++)
1645	      {
1646	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1647	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
1648	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005493 <+21651>:	movaps 0xf00(%r8),%xmm9
   0x00000000000054af <+21679>:	movaps 0x380(%rcx),%xmm5
   0x00000000000054b6 <+21686>:	mulps  %xmm5,%xmm9
   0x00000000000054c6 <+21702>:	addps  %xmm12,%xmm9
   0x00000000000054ca <+21706>:	movaps 0xf40(%r8),%xmm12
   0x00000000000054d2 <+21714>:	mulps  %xmm5,%xmm12
   0x00000000000054f9 <+21753>:	movaps 0xfc0(%r8),%xmm7
   0x0000000000005501 <+21761>:	mulps  %xmm5,%xmm7
   0x0000000000005504 <+21764>:	addps  %xmm14,%xmm12
   0x0000000000005514 <+21780>:	addps  %xmm13,%xmm7
   0x0000000000005538 <+21816>:	movaps 0xf80(%r8),%xmm14
   0x0000000000005540 <+21824>:	mulps  %xmm5,%xmm14
   0x0000000000005554 <+21844>:	addps  %xmm15,%xmm14

1649
1650	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1651	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
1652	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000545f <+21599>:	movaps 0xf10(%r8),%xmm6
   0x000000000000549f <+21663>:	movaps 0x390(%rcx),%xmm8
   0x00000000000054a7 <+21671>:	mulps  %xmm8,%xmm6
   0x00000000000054d6 <+21718>:	addps  %xmm9,%xmm6
   0x0000000000005508 <+21768>:	movaps 0xf50(%r8),%xmm14
   0x0000000000005510 <+21776>:	mulps  %xmm8,%xmm14
   0x0000000000005518 <+21784>:	movaps 0xfd0(%r8),%xmm13
   0x0000000000005520 <+21792>:	mulps  %xmm8,%xmm13
   0x0000000000005524 <+21796>:	addps  %xmm12,%xmm14
   0x0000000000005558 <+21848>:	movaps 0xf90(%r8),%xmm15
   0x0000000000005560 <+21856>:	mulps  %xmm8,%xmm15
   0x0000000000005564 <+21860>:	addps  %xmm7,%xmm13
   0x0000000000005580 <+21888>:	addps  %xmm14,%xmm15

1653
1654	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1655	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
1656	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005477 <+21623>:	movaps 0xf20(%r8),%xmm7
   0x0000000000005483 <+21635>:	movaps 0x3a0(%rcx),%xmm10
   0x000000000000548b <+21643>:	mulps  %xmm10,%xmm7
   0x00000000000054e6 <+21734>:	addps  %xmm6,%xmm7
   0x00000000000054e9 <+21737>:	movaps 0xf60(%r8),%xmm6
   0x00000000000054f1 <+21745>:	mulps  %xmm10,%xmm6
   0x0000000000005534 <+21812>:	addps  %xmm14,%xmm6
   0x0000000000005548 <+21832>:	movaps 0xfa0(%r8),%xmm6
   0x0000000000005550 <+21840>:	mulps  %xmm10,%xmm6
   0x0000000000005568 <+21864>:	movaps 0xfe0(%r8),%xmm8
   0x0000000000005570 <+21872>:	mulps  %xmm10,%xmm8
   0x0000000000005590 <+21904>:	addps  %xmm15,%xmm6
   0x00000000000055d9 <+21977>:	addps  %xmm13,%xmm8

1657
1658	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1659	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
1660	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000054be <+21694>:	movaps 0xf30(%r8),%xmm11
   0x00000000000054da <+21722>:	movaps 0x3b0(%rcx),%xmm9
   0x00000000000054e2 <+21730>:	mulps  %xmm9,%xmm11
   0x00000000000054f5 <+21749>:	addps  %xmm7,%xmm11
   0x0000000000005528 <+21800>:	movaps 0xf70(%r8),%xmm12
   0x0000000000005530 <+21808>:	mulps  %xmm9,%xmm12
   0x0000000000005544 <+21828>:	addps  %xmm6,%xmm12
   0x0000000000005574 <+21876>:	movaps 0xff0(%r8),%xmm10
   0x000000000000557c <+21884>:	mulps  %xmm9,%xmm10
   0x0000000000005584 <+21892>:	movaps 0xfb0(%r8),%xmm14
   0x000000000000558c <+21900>:	mulps  %xmm9,%xmm14
   0x00000000000055c1 <+21953>:	addps  %xmm6,%xmm14
   0x00000000000055dd <+21981>:	addps  %xmm8,%xmm10

1661	      }
1662
1663	      /* Store C(4,3) block. */
1664	      for (i = 0; i < 4; i++)
1665	      {
1666	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000005594 <+21908>:	movaps 0x20(%rsp),%xmm5
   0x0000000000005599 <+21913>:	mulps  %xmm5,%xmm11
   0x00000000000055ad <+21933>:	mulps  %xmm5,%xmm12
   0x00000000000055c5 <+21957>:	mulps  %xmm5,%xmm14
   0x00000000000055e1 <+21985>:	mulps  %xmm5,%xmm10

1667	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_43]), C_row[i]);
   0x000000000000559d <+21917>:	addps  0x380(%r9),%xmm11
   0x00000000000055b1 <+21937>:	addps  0x390(%r9),%xmm12
   0x00000000000055c9 <+21961>:	addps  0x3a0(%r9),%xmm14
   0x00000000000055e5 <+21989>:	addps  0x3b0(%r9),%xmm10

1668	        _mm_store_ps(&C[i*4+C_OFFSET_43], C_row[i]);
   0x00000000000055a5 <+21925>:	movaps %xmm11,0x380(%r9)
   0x00000000000055b9 <+21945>:	movaps %xmm12,0x390(%r9)
   0x00000000000055d1 <+21969>:	movaps %xmm14,0x3a0(%r9)
   0x00000000000055ed <+21997>:	movaps %xmm10,0x3b0(%r9)

1669	      }
1670	    }
1671
1672	    /* Reset C(4,4) matrix accumulators */
1673	    C_row[0] = _mm_setzero_ps();
1674	    C_row[1] = _mm_setzero_ps();
1675	    C_row[2] = _mm_setzero_ps();
1676	    C_row[3] = _mm_setzero_ps();
1677
1678	    if (norm_product[12][3] &&
   0x00000000000055f5 <+22005>:	comiss %xmm1,%xmm4
   0x00000000000055f8 <+22008>:	jb     0x5ab5 <stream_kernel+23221>

1679	        norm_product[13][7] &&
   0x00000000000055fe <+22014>:	comiss %xmm1,%xmm3
   0x0000000000005601 <+22017>:	jb     0x5ab5 <stream_kernel+23221>

1680	        norm_product[14][11] &&
   0x0000000000005607 <+22023>:	comiss %xmm1,%xmm2
   0x000000000000560a <+22026>:	jb     0x5ab5 <stream_kernel+23221>

1681	        norm_product[15][15])
   0x0000000000005610 <+22032>:	comiss %xmm1,%xmm0
   0x0000000000005613 <+22035>:	jb     0x5ab5 <stream_kernel+23221>

1682	    {
1683	      /* A(4,1)*B(1,4) = C(4,4). */
1684	      for (i = 0; i < 4; i++)
1685	      {
1686	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1687	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
1688	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005619 <+22041>:	movaps 0xc0(%rcx),%xmm11
   0x0000000000005639 <+22073>:	movaps 0xc00(%r8),%xmm0
   0x0000000000005659 <+22105>:	movaps 0xc40(%r8),%xmm5
   0x0000000000005679 <+22137>:	movaps 0xc80(%r8),%xmm8
   0x0000000000005691 <+22161>:	mulps  %xmm11,%xmm0
   0x000000000000569c <+22172>:	movaps 0xcc0(%r8),%xmm0
   0x00000000000056c1 <+22209>:	mulps  %xmm11,%xmm5
   0x00000000000056f1 <+22257>:	mulps  %xmm11,%xmm8
   0x0000000000005725 <+22309>:	mulps  %xmm11,%xmm0

1689
1690	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1691	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
1692	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005621 <+22049>:	movaps 0xd0(%rcx),%xmm10
   0x0000000000005641 <+22081>:	movaps 0xc10(%r8),%xmm2
   0x0000000000005661 <+22113>:	movaps 0xc50(%r8),%xmm7
   0x0000000000005681 <+22145>:	movaps 0xc90(%r8),%xmm9
   0x0000000000005695 <+22165>:	mulps  %xmm10,%xmm2
   0x0000000000005699 <+22169>:	addps  %xmm0,%xmm2
   0x00000000000056c5 <+22213>:	mulps  %xmm10,%xmm7
   0x00000000000056c9 <+22217>:	addps  %xmm5,%xmm7
   0x00000000000056f5 <+22261>:	mulps  %xmm10,%xmm9
   0x00000000000056f9 <+22265>:	addps  %xmm8,%xmm9
   0x0000000000005729 <+22313>:	movaps 0xcd0(%r8),%xmm11
   0x0000000000005731 <+22321>:	mulps  %xmm10,%xmm11
   0x000000000000573d <+22333>:	addps  %xmm0,%xmm11

1693
1694	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1695	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
1696	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005629 <+22057>:	movaps 0xe0(%rcx),%xmm13
   0x0000000000005649 <+22089>:	movaps 0xc20(%r8),%xmm4
   0x0000000000005669 <+22121>:	movaps 0xc60(%r8),%xmm6
   0x0000000000005689 <+22153>:	movaps 0xca0(%r8),%xmm15
   0x00000000000056a4 <+22180>:	mulps  %xmm13,%xmm4
   0x00000000000056a8 <+22184>:	addps  %xmm2,%xmm4
   0x00000000000056d4 <+22228>:	mulps  %xmm13,%xmm6
   0x00000000000056d8 <+22232>:	addps  %xmm7,%xmm6
   0x0000000000005705 <+22277>:	mulps  %xmm13,%xmm15
   0x0000000000005709 <+22281>:	addps  %xmm9,%xmm15
   0x0000000000005735 <+22325>:	movaps 0xce0(%r8),%xmm10
   0x0000000000005741 <+22337>:	mulps  %xmm13,%xmm10
   0x000000000000574d <+22349>:	addps  %xmm11,%xmm10

1697
1698	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1699	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
1700	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005631 <+22065>:	movaps 0xf0(%rcx),%xmm12
   0x0000000000005651 <+22097>:	movaps 0xc30(%r8),%xmm14
   0x0000000000005671 <+22129>:	movaps 0xc70(%r8),%xmm3
   0x00000000000056b2 <+22194>:	mulps  %xmm12,%xmm14
   0x00000000000056b6 <+22198>:	addps  %xmm4,%xmm14
   0x00000000000056e3 <+22243>:	mulps  %xmm12,%xmm3
   0x00000000000056e7 <+22247>:	addps  %xmm6,%xmm3
   0x00000000000056fd <+22269>:	movaps 0xcf0(%r8),%xmm8
   0x000000000000570d <+22285>:	movaps 0xcb0(%r8),%xmm9
   0x0000000000005715 <+22293>:	mulps  %xmm12,%xmm9
   0x0000000000005719 <+22297>:	addps  %xmm15,%xmm9
   0x0000000000005751 <+22353>:	mulps  %xmm12,%xmm8
   0x000000000000575d <+22365>:	addps  %xmm10,%xmm8

1701	      }
1702
1703	      /* A(4,2)*B(2,4) = C(4,4). */
1704	      for (i = 0; i < 4; i++)
1705	      {
1706	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1707	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
1708	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000056ea <+22250>:	movaps 0x1c0(%rcx),%xmm6
   0x0000000000005745 <+22341>:	movaps 0xd00(%r8),%xmm13
   0x0000000000005761 <+22369>:	mulps  %xmm6,%xmm13
   0x000000000000576c <+22380>:	addps  %xmm14,%xmm13
   0x000000000000577c <+22396>:	movaps 0xd80(%r8),%xmm11
   0x0000000000005784 <+22404>:	movaps 0xdc0(%r8),%xmm14
   0x0000000000005796 <+22422>:	movaps 0xd40(%r8),%xmm13
   0x00000000000057a2 <+22434>:	mulps  %xmm6,%xmm13
   0x00000000000057bd <+22461>:	addps  %xmm3,%xmm13
   0x00000000000057df <+22495>:	mulps  %xmm6,%xmm11
   0x00000000000057fb <+22523>:	addps  %xmm9,%xmm11
   0x000000000000581a <+22554>:	mulps  %xmm6,%xmm14
   0x000000000000585c <+22620>:	addps  %xmm8,%xmm14

1709
1710	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1711	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
1712	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000056ba <+22202>:	movaps 0x1d0(%rcx),%xmm4
   0x0000000000005755 <+22357>:	movaps 0xd10(%r8),%xmm12
   0x0000000000005770 <+22384>:	mulps  %xmm4,%xmm12
   0x000000000000578f <+22415>:	addps  %xmm13,%xmm12
   0x00000000000057c1 <+22465>:	movaps 0xd50(%r8),%xmm3
   0x00000000000057c9 <+22473>:	mulps  %xmm4,%xmm3
   0x00000000000057cc <+22476>:	addps  %xmm13,%xmm3
   0x0000000000005802 <+22530>:	movaps 0xd90(%r8),%xmm9
   0x000000000000580a <+22538>:	mulps  %xmm4,%xmm9
   0x000000000000580e <+22542>:	addps  %xmm11,%xmm9
   0x000000000000581e <+22558>:	movaps 0xdd0(%r8),%xmm6
   0x000000000000582a <+22570>:	mulps  %xmm4,%xmm6
   0x000000000000586c <+22636>:	addps  %xmm14,%xmm6

1713
1714	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1715	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
1716	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000056db <+22235>:	movaps 0xd20(%r8),%xmm7
   0x000000000000571d <+22301>:	movaps 0xda0(%r8),%xmm15
   0x0000000000005765 <+22373>:	movaps 0x1e0(%rcx),%xmm0
   0x0000000000005774 <+22388>:	movaps 0xd60(%r8),%xmm10
   0x000000000000578c <+22412>:	mulps  %xmm0,%xmm7
   0x000000000000579e <+22430>:	addps  %xmm12,%xmm7
   0x00000000000057b1 <+22449>:	mulps  %xmm0,%xmm10
   0x00000000000057db <+22491>:	addps  %xmm3,%xmm10
   0x00000000000057ef <+22511>:	mulps  %xmm0,%xmm15
   0x0000000000005826 <+22566>:	addps  %xmm9,%xmm15
   0x000000000000582d <+22573>:	movaps 0xde0(%r8),%xmm4
   0x0000000000005841 <+22593>:	mulps  %xmm0,%xmm4
   0x000000000000587c <+22652>:	addps  %xmm6,%xmm4

1717
1718	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1719	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
1720	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000056ab <+22187>:	movaps 0x1f0(%rcx),%xmm2
   0x00000000000056cc <+22220>:	movaps 0xd30(%r8),%xmm5
   0x0000000000005793 <+22419>:	mulps  %xmm2,%xmm5
   0x00000000000057ae <+22446>:	addps  %xmm7,%xmm5
   0x00000000000057b5 <+22453>:	movaps 0xd70(%r8),%xmm7
   0x00000000000057d0 <+22480>:	mulps  %xmm2,%xmm7
   0x00000000000057e3 <+22499>:	movaps 0xdb0(%r8),%xmm3
   0x00000000000057eb <+22507>:	addps  %xmm10,%xmm7
   0x00000000000057ff <+22527>:	mulps  %xmm2,%xmm3
   0x0000000000005835 <+22581>:	addps  %xmm15,%xmm3
   0x0000000000005839 <+22585>:	movaps 0xdf0(%r8),%xmm15
   0x000000000000584c <+22604>:	mulps  %xmm2,%xmm15
   0x0000000000005897 <+22679>:	addps  %xmm4,%xmm15

1721	      }
1722
1723	      /* A(4,3)*B(3,4) = C(4,4). */
1724	      for (i = 0; i < 4; i++)
1725	      {
1726	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1727	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
1728	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005844 <+22596>:	movaps 0xe00(%r8),%xmm0
   0x0000000000005870 <+22640>:	movaps 0x2c0(%rcx),%xmm14
   0x0000000000005878 <+22648>:	mulps  %xmm14,%xmm0
   0x00000000000058a7 <+22695>:	addps  %xmm5,%xmm0
   0x00000000000058b9 <+22713>:	movaps 0xe40(%r8),%xmm0
   0x00000000000058c1 <+22721>:	mulps  %xmm14,%xmm0
   0x00000000000058e3 <+22755>:	addps  %xmm7,%xmm0
   0x00000000000058e6 <+22758>:	movaps 0xe80(%r8),%xmm7
   0x00000000000058ee <+22766>:	mulps  %xmm14,%xmm7
   0x0000000000005915 <+22805>:	addps  %xmm3,%xmm7
   0x0000000000005942 <+22850>:	movaps 0xec0(%r8),%xmm3
   0x000000000000594a <+22858>:	mulps  %xmm14,%xmm3
   0x000000000000597d <+22909>:	addps  %xmm15,%xmm3

1729
1730	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1731	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
1732	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005812 <+22546>:	movaps 0x2d0(%rcx),%xmm11
   0x0000000000005850 <+22608>:	movaps 0xe10(%r8),%xmm2
   0x0000000000005858 <+22616>:	mulps  %xmm11,%xmm2
   0x00000000000058b6 <+22710>:	addps  %xmm0,%xmm2
   0x00000000000058c8 <+22728>:	movaps 0xe50(%r8),%xmm2
   0x00000000000058d0 <+22736>:	mulps  %xmm11,%xmm2
   0x00000000000058f2 <+22770>:	addps  %xmm0,%xmm2
   0x0000000000005918 <+22808>:	movaps 0xe90(%r8),%xmm3
   0x0000000000005920 <+22816>:	mulps  %xmm11,%xmm3
   0x0000000000005924 <+22820>:	addps  %xmm7,%xmm3
   0x000000000000594e <+22862>:	movaps 0xed0(%r8),%xmm14
   0x0000000000005956 <+22870>:	mulps  %xmm11,%xmm14
   0x000000000000598d <+22925>:	addps  %xmm3,%xmm14

1733
1734	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1735	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
1736	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005864 <+22628>:	movaps 0x2e0(%rcx),%xmm9
   0x0000000000005883 <+22659>:	movaps 0xea0(%r8),%xmm6
   0x0000000000005893 <+22675>:	mulps  %xmm9,%xmm6
   0x000000000000589f <+22687>:	movaps 0xe60(%r8),%xmm4
   0x00000000000058aa <+22698>:	movaps 0xe20(%r8),%xmm5
   0x00000000000058b2 <+22706>:	mulps  %xmm9,%xmm5
   0x00000000000058c5 <+22725>:	addps  %xmm2,%xmm5
   0x00000000000058d8 <+22744>:	mulps  %xmm9,%xmm4
   0x00000000000058fd <+22781>:	addps  %xmm2,%xmm4
   0x000000000000593f <+22847>:	addps  %xmm3,%xmm6
   0x000000000000595a <+22874>:	movaps 0xee0(%r8),%xmm11
   0x0000000000005962 <+22882>:	mulps  %xmm9,%xmm11
   0x000000000000599c <+22940>:	addps  %xmm14,%xmm11

1737
1738	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1739	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
1740	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000057a6 <+22438>:	movaps 0xe30(%r8),%xmm12
   0x00000000000057d3 <+22483>:	movaps 0xe70(%r8),%xmm13
   0x00000000000057f3 <+22515>:	movaps 0x2f0(%rcx),%xmm10
   0x0000000000005860 <+22624>:	mulps  %xmm10,%xmm12
   0x000000000000587f <+22655>:	mulps  %xmm10,%xmm13
   0x000000000000588b <+22667>:	movaps 0xeb0(%r8),%xmm8
   0x000000000000589b <+22683>:	mulps  %xmm10,%xmm8
   0x00000000000058d4 <+22740>:	addps  %xmm5,%xmm12
   0x000000000000590a <+22794>:	addps  %xmm4,%xmm13
   0x0000000000005927 <+22823>:	movaps 0xef0(%r8),%xmm7
   0x000000000000592f <+22831>:	mulps  %xmm10,%xmm7
   0x000000000000596e <+22894>:	addps  %xmm6,%xmm8
   0x00000000000059a0 <+22944>:	addps  %xmm11,%xmm7

1741	      }
1742
1743	      /* A(4,4)*B(4,4) = C(4,4). */
1744	      for (i = 0; i < 4; i++)
1745	      {
1746	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1747	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
1748	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005966 <+22886>:	movaps 0xf00(%r8),%xmm9
   0x0000000000005972 <+22898>:	movaps 0x3c0(%rcx),%xmm6
   0x0000000000005979 <+22905>:	mulps  %xmm6,%xmm9
   0x0000000000005981 <+22913>:	movaps 0xf40(%r8),%xmm15
   0x0000000000005989 <+22921>:	mulps  %xmm6,%xmm15
   0x00000000000059b0 <+22960>:	addps  %xmm12,%xmm9
   0x00000000000059d4 <+22996>:	movaps 0xf80(%r8),%xmm10
   0x00000000000059dc <+23004>:	mulps  %xmm6,%xmm10
   0x00000000000059e4 <+23012>:	addps  %xmm13,%xmm15
   0x00000000000059f4 <+23028>:	addps  %xmm8,%xmm10
   0x0000000000005a20 <+23072>:	movaps 0xfc0(%r8),%xmm10
   0x0000000000005a28 <+23080>:	mulps  %xmm6,%xmm10
   0x0000000000005a94 <+23188>:	addps  %xmm7,%xmm10

1749
1750	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1751	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
1752	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000058dc <+22748>:	movaps 0x3d0(%rcx),%xmm5
   0x0000000000005933 <+22835>:	movaps 0xf10(%r8),%xmm10
   0x000000000000593b <+22843>:	mulps  %xmm5,%xmm10
   0x00000000000059c0 <+22976>:	addps  %xmm9,%xmm10
   0x00000000000059e8 <+23016>:	movaps 0xf50(%r8),%xmm13
   0x00000000000059f0 <+23024>:	mulps  %xmm5,%xmm13
   0x00000000000059f8 <+23032>:	movaps 0xf90(%r8),%xmm8
   0x0000000000005a00 <+23040>:	mulps  %xmm5,%xmm8
   0x0000000000005a04 <+23044>:	addps  %xmm15,%xmm13
   0x0000000000005a1c <+23068>:	addps  %xmm10,%xmm8
   0x0000000000005a2c <+23084>:	movaps 0xfd0(%r8),%xmm6
   0x0000000000005a34 <+23092>:	mulps  %xmm5,%xmm6
   0x0000000000005a98 <+23192>:	addps  %xmm10,%xmm6

1753
1754	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1755	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
1756	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000590e <+22798>:	movaps 0x3e0(%rcx),%xmm4
   0x00000000000059b4 <+22964>:	movaps 0xf20(%r8),%xmm12
   0x00000000000059bc <+22972>:	mulps  %xmm4,%xmm12
   0x00000000000059c4 <+22980>:	movaps 0xf60(%r8),%xmm9
   0x00000000000059cc <+22988>:	mulps  %xmm4,%xmm9
   0x00000000000059d0 <+22992>:	addps  %xmm10,%xmm12
   0x0000000000005a08 <+23048>:	addps  %xmm13,%xmm9
   0x0000000000005a10 <+23056>:	movaps 0xfa0(%r8),%xmm9
   0x0000000000005a18 <+23064>:	mulps  %xmm4,%xmm9
   0x0000000000005a37 <+23095>:	movaps 0xfe0(%r8),%xmm5
   0x0000000000005a3f <+23103>:	mulps  %xmm4,%xmm5
   0x0000000000005a4d <+23117>:	addps  %xmm8,%xmm9
   0x0000000000005a9c <+23196>:	addps  %xmm6,%xmm5

1757
1758	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1759	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
1760	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000058f5 <+22773>:	movaps 0xf70(%r8),%xmm0
   0x0000000000005900 <+22784>:	movaps 0x3f0(%rcx),%xmm2
   0x0000000000005907 <+22791>:	mulps  %xmm2,%xmm0
   0x0000000000005991 <+22929>:	movaps 0xf30(%r8),%xmm3
   0x0000000000005999 <+22937>:	mulps  %xmm2,%xmm3
   0x00000000000059a4 <+22948>:	movaps 0xfb0(%r8),%xmm11
   0x00000000000059ac <+22956>:	mulps  %xmm2,%xmm11
   0x00000000000059e0 <+23008>:	addps  %xmm12,%xmm3
   0x0000000000005a0c <+23052>:	addps  %xmm9,%xmm0
   0x0000000000005a42 <+23106>:	movaps 0xff0(%r8),%xmm4
   0x0000000000005a4a <+23114>:	mulps  %xmm2,%xmm4
   0x0000000000005a7c <+23164>:	addps  %xmm9,%xmm11
   0x0000000000005a9f <+23199>:	addps  %xmm5,%xmm4

1761	      }
1762
1763	      /* Store C(4,4) block. */
1764	      for (i = 0; i < 4; i++)
1765	      {
1766	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000005a51 <+23121>:	movaps 0x20(%rsp),%xmm2
   0x0000000000005a56 <+23126>:	mulps  %xmm2,%xmm3
   0x0000000000005a69 <+23145>:	mulps  %xmm2,%xmm0
   0x0000000000005a80 <+23168>:	mulps  %xmm2,%xmm11
   0x0000000000005aa2 <+23202>:	mulps  %xmm2,%xmm4

1767	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_44]), C_row[i]);
   0x0000000000005a59 <+23129>:	addps  0x3c0(%r9),%xmm3
   0x0000000000005a6c <+23148>:	addps  0x3d0(%r9),%xmm0
   0x0000000000005a84 <+23172>:	addps  0x3e0(%r9),%xmm11
   0x0000000000005aa5 <+23205>:	addps  0x3f0(%r9),%xmm4

1768	        _mm_store_ps(&C[i*4+C_OFFSET_44], C_row[i]);
   0x0000000000005a61 <+23137>:	movaps %xmm3,0x3c0(%r9)
   0x0000000000005a74 <+23156>:	movaps %xmm0,0x3d0(%r9)
   0x0000000000005a8c <+23180>:	movaps %xmm11,0x3e0(%r9)
   0x0000000000005aad <+23213>:	movaps %xmm4,0x3f0(%r9)

1769	      }
1770	    }
1771	  }
1772	}
   0x0000000000005ac3 <+23235>:	add    $0x208,%rsp
   0x0000000000005aca <+23242>:	retq
   0x0000000000005acb <+23243>:	nopl   0x0(%rax,%rax,1)

End of assembler dump.
