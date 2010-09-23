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
91	  /* Divide number of stream elements by 64 to simulate stride of 64. */
92	  max_stream_index = number_stream_elements/64;
   0x0000000000000009 <+9>:	shr    $0x6,%edi

93
94	  alpha_row = _mm_set1_ps(alpha);
   0x0000000000000000 <+0>:	movaps %xmm0,%xmm2
   0x0000000000000003 <+3>:	shufps $0x0,%xmm2,%xmm2

95
96	  for (stream_index = 0; stream_index < max_stream_index; stream_index++)
   0x0000000000000007 <+7>:	xor    %edx,%edx
   0x000000000000000c <+12>:	test   %edi,%edi
   0x000000000000000e <+14>:	jbe    0x500d <stream_kernel+20493>
   0x0000000000000014 <+20>:	xor    %eax,%eax
   0x0000000000004fff <+20479>:	inc    %eax
   0x0000000000005001 <+20481>:	mov    %eax,%edx
   0x0000000000005003 <+20483>:	mov    %edx,%eax
   0x0000000000005005 <+20485>:	cmp    %edi,%eax
   0x0000000000005007 <+20487>:	jb     0x16 <stream_kernel+22>

97	  {
98	    /* Load pointers to matrix data blocks. */
99	    A = multiply_stream[stream_index].A_block;
   0x0000000000000016 <+22>:	imul   $0x98,%rdx,%rdx
   0x000000000000001d <+29>:	mov    (%rsi,%rdx,1),%r8

100	    B = multiply_stream[stream_index].B_block;
   0x0000000000000021 <+33>:	mov    0x8(%rsi,%rdx,1),%rcx

101	    C = multiply_stream[stream_index].C_block;
   0x0000000000000026 <+38>:	mov    0x10(%rsi,%rdx,1),%r9

102	    norm = multiply_stream[stream_index].norm;
103
104	    /* Reset C(1,1) matrix accumulators */
105	    C_row[0] = _mm_setzero_ps();
106	    C_row[1] = _mm_setzero_ps();
107	    C_row[2] = _mm_setzero_ps();
108	    C_row[3] = _mm_setzero_ps();
109
110	    if (norm[0]*norm[16] >= tolerance &&
   0x000000000000002b <+43>:	movss  0x18(%rdx,%rsi,1),%xmm7
   0x0000000000000031 <+49>:	movss  0x58(%rdx,%rsi,1),%xmm6
   0x0000000000000037 <+55>:	movaps %xmm7,%xmm0
   0x000000000000003a <+58>:	mulss  %xmm6,%xmm0
   0x000000000000003e <+62>:	comiss %xmm1,%xmm0
   0x0000000000000041 <+65>:	jb     0x4f7 <stream_kernel+1271>

111	        norm[1]*norm[20] >= tolerance &&
   0x0000000000000047 <+71>:	movss  0x1c(%rdx,%rsi,1),%xmm0
   0x000000000000004d <+77>:	mulss  0x68(%rdx,%rsi,1),%xmm0
   0x0000000000000053 <+83>:	comiss %xmm1,%xmm0
   0x0000000000000056 <+86>:	jb     0x4f7 <stream_kernel+1271>

112	        norm[2]*norm[24] >= tolerance &&
   0x000000000000005c <+92>:	movss  0x20(%rdx,%rsi,1),%xmm0
   0x0000000000000062 <+98>:	mulss  0x78(%rdx,%rsi,1),%xmm0
   0x0000000000000068 <+104>:	comiss %xmm1,%xmm0
   0x000000000000006b <+107>:	jb     0x4f7 <stream_kernel+1271>

113	        norm[3]*norm[28] >= tolerance)
   0x0000000000000071 <+113>:	movss  0x24(%rdx,%rsi,1),%xmm0
   0x0000000000000077 <+119>:	mulss  0x88(%rdx,%rsi,1),%xmm0
   0x0000000000000080 <+128>:	comiss %xmm1,%xmm0
   0x0000000000000083 <+131>:	jb     0x4f7 <stream_kernel+1271>

114	    {
115	      /* A(1,1)*B(1,1) = C(1,1). */
116	      for (i = 0; i < 4; i++)
117	      {
118	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
119	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
   0x0000000000000089 <+137>:	movaps (%rcx),%xmm4

120	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000009b <+155>:	movaps (%r8),%xmm0
   0x00000000000000ae <+174>:	movaps 0x40(%r8),%xmm11
   0x00000000000000bd <+189>:	mulps  %xmm4,%xmm0
   0x00000000000000c4 <+196>:	movaps 0xc0(%r8),%xmm8
   0x00000000000000dd <+221>:	movaps 0x80(%r8),%xmm12
   0x00000000000000f5 <+245>:	mulps  %xmm4,%xmm11
   0x0000000000000128 <+296>:	mulps  %xmm4,%xmm12
   0x000000000000015c <+348>:	mulps  %xmm4,%xmm8

121
122	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
123	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
124	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000008c <+140>:	movaps 0x10(%rcx),%xmm15
   0x000000000000009f <+159>:	movaps 0x10(%r8),%xmm12
   0x00000000000000b3 <+179>:	movaps 0x50(%r8),%xmm9
   0x00000000000000c0 <+192>:	mulps  %xmm15,%xmm12
   0x00000000000000cc <+204>:	addps  %xmm0,%xmm12
   0x00000000000000ed <+237>:	movaps 0x90(%r8),%xmm13
   0x00000000000000f9 <+249>:	mulps  %xmm15,%xmm9
   0x00000000000000fd <+253>:	addps  %xmm11,%xmm9
   0x000000000000012c <+300>:	mulps  %xmm15,%xmm13
   0x0000000000000130 <+304>:	addps  %xmm12,%xmm13
   0x0000000000000160 <+352>:	movaps 0xd0(%r8),%xmm4
   0x0000000000000168 <+360>:	mulps  %xmm15,%xmm4
   0x0000000000000174 <+372>:	addps  %xmm8,%xmm4

125
126	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
127	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
128	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000091 <+145>:	movaps 0x20(%rcx),%xmm10
   0x00000000000000a4 <+164>:	movaps 0x20(%r8),%xmm13
   0x00000000000000b8 <+184>:	movaps 0x60(%r8),%xmm3
   0x00000000000000d0 <+208>:	mulps  %xmm10,%xmm13
   0x00000000000000d9 <+217>:	addps  %xmm12,%xmm13
   0x0000000000000101 <+257>:	movaps 0xa0(%r8),%xmm11
   0x0000000000000109 <+265>:	mulps  %xmm10,%xmm3
   0x000000000000010d <+269>:	addps  %xmm9,%xmm3
   0x000000000000013c <+316>:	mulps  %xmm10,%xmm11
   0x0000000000000140 <+320>:	addps  %xmm13,%xmm11
   0x000000000000016c <+364>:	movaps 0xe0(%r8),%xmm15
   0x0000000000000178 <+376>:	mulps  %xmm10,%xmm15
   0x0000000000000184 <+388>:	addps  %xmm4,%xmm15

129
130	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
131	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
132	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000096 <+150>:	movaps 0x30(%rcx),%xmm14
   0x00000000000000a9 <+169>:	movaps 0x30(%r8),%xmm5
   0x00000000000000d4 <+212>:	movaps 0x70(%r8),%xmm0
   0x00000000000000e5 <+229>:	mulps  %xmm14,%xmm5
   0x00000000000000e9 <+233>:	addps  %xmm13,%xmm5
   0x0000000000000119 <+281>:	mulps  %xmm14,%xmm0
   0x000000000000011d <+285>:	addps  %xmm3,%xmm0
   0x0000000000000120 <+288>:	movaps 0xb0(%r8),%xmm3
   0x000000000000014c <+332>:	mulps  %xmm14,%xmm3
   0x0000000000000150 <+336>:	addps  %xmm11,%xmm3
   0x0000000000000154 <+340>:	movaps 0xf0(%r8),%xmm11
   0x0000000000000190 <+400>:	mulps  %xmm14,%xmm11
   0x00000000000001a0 <+416>:	addps  %xmm15,%xmm11

133	      }
134
135	      /* A(1,2)*B(2,1) = C(1,1). */
136	      for (i = 0; i < 4; i++)
137	      {
138	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
139	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
140	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000144 <+324>:	movaps 0x100(%rcx),%xmm13
   0x000000000000017c <+380>:	movaps 0x100(%r8),%xmm10
   0x0000000000000188 <+392>:	movaps 0x140(%r8),%xmm4
   0x000000000000019c <+412>:	mulps  %xmm13,%xmm10
   0x00000000000001bc <+444>:	addps  %xmm5,%xmm10
   0x00000000000001c0 <+448>:	mulps  %xmm13,%xmm4
   0x00000000000001dc <+476>:	addps  %xmm0,%xmm4
   0x00000000000001ef <+495>:	movaps 0x180(%r8),%xmm14
   0x00000000000001f7 <+503>:	mulps  %xmm13,%xmm14
   0x000000000000022e <+558>:	movaps 0x1c0(%r8),%xmm10
   0x000000000000023a <+570>:	addps  %xmm3,%xmm14
   0x000000000000024e <+590>:	mulps  %xmm13,%xmm10
   0x0000000000000288 <+648>:	addps  %xmm11,%xmm10

141
142	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
143	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
144	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000194 <+404>:	movaps 0x110(%r8),%xmm14
   0x00000000000001a8 <+424>:	movaps 0x110(%rcx),%xmm8
   0x00000000000001b8 <+440>:	mulps  %xmm8,%xmm14
   0x00000000000001c4 <+452>:	movaps 0x150(%r8),%xmm5
   0x00000000000001cc <+460>:	addps  %xmm10,%xmm14
   0x00000000000001df <+479>:	mulps  %xmm8,%xmm5
   0x000000000000020b <+523>:	addps  %xmm4,%xmm5
   0x000000000000023e <+574>:	movaps 0x190(%r8),%xmm3
   0x0000000000000246 <+582>:	mulps  %xmm8,%xmm3
   0x000000000000024a <+586>:	addps  %xmm14,%xmm3
   0x0000000000000252 <+594>:	movaps 0x1d0(%r8),%xmm13
   0x0000000000000265 <+613>:	mulps  %xmm8,%xmm13
   0x00000000000002a4 <+676>:	addps  %xmm10,%xmm13

145
146	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
147	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
148	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000001b0 <+432>:	movaps 0x120(%rcx),%xmm15
   0x00000000000001d0 <+464>:	movaps 0x120(%r8),%xmm10
   0x00000000000001d8 <+472>:	mulps  %xmm15,%xmm10
   0x00000000000001e3 <+483>:	movaps 0x1a0(%r8),%xmm0
   0x00000000000001eb <+491>:	addps  %xmm14,%xmm10
   0x00000000000001ff <+511>:	movaps 0x160(%r8),%xmm10
   0x0000000000000207 <+519>:	mulps  %xmm15,%xmm10
   0x000000000000021a <+538>:	addps  %xmm5,%xmm10
   0x0000000000000226 <+550>:	mulps  %xmm15,%xmm0
   0x000000000000025a <+602>:	addps  %xmm3,%xmm0
   0x000000000000028c <+652>:	movaps 0x1e0(%r8),%xmm11
   0x0000000000000294 <+660>:	mulps  %xmm15,%xmm11
   0x00000000000002b8 <+696>:	addps  %xmm13,%xmm11

149
150	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
151	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
152	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000111 <+273>:	movaps 0x130(%rcx),%xmm9
   0x0000000000000134 <+308>:	movaps 0x130(%r8),%xmm12
   0x00000000000001a4 <+420>:	mulps  %xmm9,%xmm12
   0x00000000000001fb <+507>:	addps  %xmm10,%xmm12
   0x000000000000020e <+526>:	movaps 0x170(%r8),%xmm4
   0x0000000000000216 <+534>:	mulps  %xmm9,%xmm4
   0x000000000000021e <+542>:	movaps 0x1b0(%r8),%xmm5
   0x000000000000022a <+554>:	addps  %xmm10,%xmm4
   0x0000000000000236 <+566>:	mulps  %xmm9,%xmm5
   0x000000000000025d <+605>:	movaps 0x1f0(%r8),%xmm3
   0x0000000000000271 <+625>:	mulps  %xmm9,%xmm3
   0x000000000000027d <+637>:	addps  %xmm0,%xmm5
   0x00000000000002c8 <+712>:	addps  %xmm11,%xmm3

153	      }
154
155	      /* A(1,3)*B(3,1) = C(1,1). */
156	      for (i = 0; i < 4; i++)
157	      {
158	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
159	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
160	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000280 <+640>:	movaps 0x200(%rcx),%xmm14
   0x0000000000000298 <+664>:	movaps 0x200(%r8),%xmm15
   0x00000000000002a0 <+672>:	mulps  %xmm14,%xmm15
   0x00000000000002d8 <+728>:	addps  %xmm12,%xmm15
   0x00000000000002f0 <+752>:	movaps 0x240(%r8),%xmm15
   0x00000000000002fc <+764>:	mulps  %xmm14,%xmm15
   0x0000000000000318 <+792>:	addps  %xmm4,%xmm15
   0x000000000000032c <+812>:	movaps 0x280(%r8),%xmm15
   0x0000000000000334 <+820>:	mulps  %xmm14,%xmm15
   0x000000000000033b <+827>:	movaps 0x2c0(%r8),%xmm4
   0x000000000000034b <+843>:	mulps  %xmm14,%xmm4
   0x0000000000000357 <+855>:	addps  %xmm5,%xmm15
   0x0000000000000367 <+871>:	addps  %xmm3,%xmm4

161
162	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
163	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
164	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000275 <+629>:	movaps 0x210(%r8),%xmm9
   0x00000000000002cc <+716>:	movaps 0x210(%rcx),%xmm11
   0x00000000000002d4 <+724>:	mulps  %xmm11,%xmm9
   0x00000000000002e8 <+744>:	addps  %xmm15,%xmm9
   0x000000000000031c <+796>:	movaps 0x250(%r8),%xmm4
   0x0000000000000324 <+804>:	mulps  %xmm11,%xmm4
   0x0000000000000328 <+808>:	addps  %xmm15,%xmm4
   0x000000000000034f <+847>:	movaps 0x2d0(%r8),%xmm14
   0x000000000000035b <+859>:	movaps 0x290(%r8),%xmm5
   0x0000000000000363 <+867>:	mulps  %xmm11,%xmm5
   0x000000000000036a <+874>:	mulps  %xmm11,%xmm14
   0x0000000000000376 <+886>:	addps  %xmm15,%xmm5
   0x0000000000000386 <+902>:	addps  %xmm4,%xmm14

165
166	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
167	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
168	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000269 <+617>:	movaps 0x220(%r8),%xmm8
   0x00000000000002b0 <+688>:	movaps 0x260(%r8),%xmm0
   0x00000000000002dc <+732>:	movaps 0x220(%rcx),%xmm12
   0x00000000000002e4 <+740>:	mulps  %xmm12,%xmm8
   0x00000000000002ec <+748>:	mulps  %xmm12,%xmm0
   0x00000000000002f8 <+760>:	addps  %xmm9,%xmm8
   0x0000000000000338 <+824>:	addps  %xmm4,%xmm0
   0x000000000000036e <+878>:	movaps 0x2e0(%r8),%xmm11
   0x000000000000037a <+890>:	movaps 0x2a0(%r8),%xmm15
   0x0000000000000382 <+898>:	mulps  %xmm12,%xmm15
   0x000000000000038a <+906>:	mulps  %xmm12,%xmm11
   0x00000000000003b2 <+946>:	addps  %xmm5,%xmm15
   0x00000000000003d0 <+976>:	addps  %xmm14,%xmm11

169
170	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
171	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
172	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000002a8 <+680>:	movaps 0x230(%r8),%xmm10
   0x00000000000002bc <+700>:	movaps 0x230(%rcx),%xmm13
   0x00000000000002c4 <+708>:	mulps  %xmm13,%xmm10
   0x0000000000000300 <+768>:	movaps 0x2b0(%r8),%xmm9
   0x0000000000000308 <+776>:	addps  %xmm8,%xmm10
   0x000000000000030c <+780>:	movaps 0x270(%r8),%xmm8
   0x0000000000000314 <+788>:	mulps  %xmm13,%xmm8
   0x0000000000000343 <+835>:	mulps  %xmm13,%xmm9
   0x0000000000000347 <+839>:	addps  %xmm0,%xmm8
   0x000000000000038e <+910>:	movaps 0x2f0(%r8),%xmm4
   0x00000000000003ae <+942>:	mulps  %xmm13,%xmm4
   0x00000000000003c1 <+961>:	addps  %xmm15,%xmm9
   0x00000000000003df <+991>:	addps  %xmm11,%xmm4

173	      }
174
175	      /* A(1,4)*B(4,1) = C(1,1). */
176	      for (i = 0; i < 4; i++)
177	      {
178	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
179	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
180	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000396 <+918>:	movaps 0x300(%r8),%xmm12
   0x000000000000039e <+926>:	movaps 0x340(%r8),%xmm3
   0x00000000000003a6 <+934>:	movaps 0x380(%r8),%xmm0
   0x00000000000003b6 <+950>:	movaps 0x300(%rcx),%xmm5
   0x00000000000003bd <+957>:	mulps  %xmm5,%xmm12
   0x00000000000003cd <+973>:	mulps  %xmm5,%xmm3
   0x00000000000003dc <+988>:	mulps  %xmm5,%xmm0
   0x00000000000003f7 <+1015>:	addps  %xmm10,%xmm12
   0x0000000000000427 <+1063>:	addps  %xmm8,%xmm3
   0x0000000000000473 <+1139>:	addps  %xmm9,%xmm0
   0x0000000000000497 <+1175>:	movaps 0x3c0(%r8),%xmm0
   0x000000000000049f <+1183>:	mulps  %xmm5,%xmm0
   0x00000000000004d2 <+1234>:	addps  %xmm4,%xmm0

181
182	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
183	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
184	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000003c5 <+965>:	movaps 0x310(%rcx),%xmm15
   0x00000000000003e3 <+995>:	movaps 0x310(%r8),%xmm13
   0x00000000000003eb <+1003>:	movaps 0x350(%r8),%xmm11
   0x00000000000003f3 <+1011>:	mulps  %xmm15,%xmm13
   0x0000000000000407 <+1031>:	addps  %xmm12,%xmm13
   0x000000000000040b <+1035>:	mulps  %xmm15,%xmm11
   0x0000000000000447 <+1095>:	addps  %xmm3,%xmm11
   0x0000000000000453 <+1107>:	movaps 0x3d0(%r8),%xmm3
   0x000000000000045f <+1119>:	mulps  %xmm15,%xmm3
   0x0000000000000477 <+1143>:	movaps 0x390(%r8),%xmm9
   0x000000000000047f <+1151>:	mulps  %xmm15,%xmm9
   0x0000000000000493 <+1171>:	addps  %xmm0,%xmm9
   0x00000000000004e1 <+1249>:	addps  %xmm0,%xmm3

185
186	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
187	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
188	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000003d4 <+980>:	movaps 0x320(%rcx),%xmm14
   0x00000000000003fb <+1019>:	movaps 0x320(%r8),%xmm10
   0x0000000000000403 <+1027>:	mulps  %xmm14,%xmm10
   0x0000000000000417 <+1047>:	addps  %xmm13,%xmm10
   0x000000000000042b <+1067>:	movaps 0x3a0(%r8),%xmm8
   0x0000000000000433 <+1075>:	mulps  %xmm14,%xmm8
   0x000000000000043b <+1083>:	movaps 0x360(%r8),%xmm10
   0x0000000000000443 <+1091>:	mulps  %xmm14,%xmm10
   0x0000000000000463 <+1123>:	addps  %xmm11,%xmm10
   0x00000000000004aa <+1194>:	addps  %xmm9,%xmm8
   0x00000000000004d5 <+1237>:	movaps 0x3e0(%r8),%xmm4
   0x00000000000004dd <+1245>:	mulps  %xmm14,%xmm4
   0x00000000000004e4 <+1252>:	addps  %xmm3,%xmm4

189
190	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
191	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
192	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000040f <+1039>:	movaps 0x330(%rcx),%xmm12
   0x000000000000041b <+1051>:	movaps 0x330(%r8),%xmm13
   0x0000000000000423 <+1059>:	mulps  %xmm12,%xmm13
   0x0000000000000437 <+1079>:	addps  %xmm10,%xmm13
   0x0000000000000467 <+1127>:	movaps 0x370(%r8),%xmm11
   0x000000000000046f <+1135>:	mulps  %xmm12,%xmm11
   0x0000000000000483 <+1155>:	addps  %xmm10,%xmm11
   0x0000000000000487 <+1159>:	movaps 0x3b0(%r8),%xmm10
   0x000000000000048f <+1167>:	mulps  %xmm12,%xmm10
   0x00000000000004a2 <+1186>:	movaps 0x3f0(%r8),%xmm5
   0x00000000000004ae <+1198>:	mulps  %xmm12,%xmm5
   0x00000000000004b2 <+1202>:	addps  %xmm8,%xmm10
   0x00000000000004e7 <+1255>:	addps  %xmm4,%xmm5

193	      }
194
195	      /* Store C(1,1) block. */
196	      for (i = 0; i < 4; i++)
197	      {
198	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x000000000000044b <+1099>:	mulps  %xmm2,%xmm13
   0x00000000000004b6 <+1206>:	mulps  %xmm2,%xmm11
   0x00000000000004c4 <+1220>:	mulps  %xmm2,%xmm10
   0x00000000000004ea <+1258>:	mulps  %xmm2,%xmm5

199	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_11]), C_row[i]);
   0x000000000000044f <+1103>:	addps  (%r9),%xmm13
   0x00000000000004ba <+1210>:	addps  0x10(%r9),%xmm11
   0x00000000000004c8 <+1224>:	addps  0x20(%r9),%xmm10
   0x00000000000004ed <+1261>:	addps  0x30(%r9),%xmm5

200	        _mm_store_ps(&C[i*4+C_OFFSET_11], C_row[i]);
   0x000000000000045b <+1115>:	movaps %xmm13,(%r9)
   0x00000000000004bf <+1215>:	movaps %xmm11,0x10(%r9)
   0x00000000000004cd <+1229>:	movaps %xmm10,0x20(%r9)
   0x00000000000004f2 <+1266>:	movaps %xmm5,0x30(%r9)

201	      }
202	    }
203
204	    /* Reset C(1,2) matrix accumulators */
205	    C_row[0] = _mm_setzero_ps();
206	    C_row[1] = _mm_setzero_ps();
207	    C_row[2] = _mm_setzero_ps();
208	    C_row[3] = _mm_setzero_ps();
209
210	    if (norm[0]*norm[17] >= tolerance &&
   0x00000000000004f7 <+1271>:	movss  0x5c(%rdx,%rsi,1),%xmm5
   0x00000000000004fd <+1277>:	movaps %xmm7,%xmm0
   0x0000000000000500 <+1280>:	mulss  %xmm5,%xmm0
   0x0000000000000504 <+1284>:	comiss %xmm1,%xmm0
   0x0000000000000507 <+1287>:	jb     0x9c7 <stream_kernel+2503>

211	        norm[1]*norm[21] >= tolerance &&
   0x000000000000050d <+1293>:	movss  0x1c(%rdx,%rsi,1),%xmm0
   0x0000000000000513 <+1299>:	mulss  0x6c(%rdx,%rsi,1),%xmm0
   0x0000000000000519 <+1305>:	comiss %xmm1,%xmm0
   0x000000000000051c <+1308>:	jb     0x9c7 <stream_kernel+2503>

212	        norm[2]*norm[25] >= tolerance &&
   0x0000000000000522 <+1314>:	movss  0x20(%rdx,%rsi,1),%xmm0
   0x0000000000000528 <+1320>:	mulss  0x7c(%rdx,%rsi,1),%xmm0
   0x000000000000052e <+1326>:	comiss %xmm1,%xmm0
   0x0000000000000531 <+1329>:	jb     0x9c7 <stream_kernel+2503>

213	        norm[3]*norm[29] >= tolerance)
   0x0000000000000537 <+1335>:	movss  0x24(%rdx,%rsi,1),%xmm0
   0x000000000000053d <+1341>:	mulss  0x8c(%rdx,%rsi,1),%xmm0
   0x0000000000000546 <+1350>:	comiss %xmm1,%xmm0
   0x0000000000000549 <+1353>:	jb     0x9c7 <stream_kernel+2503>

214	    {
215	      /* A(1,1)*B(1,2) = C(1,2). */
216	      for (i = 0; i < 4; i++)
217	      {
218	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
219	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
220	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000054f <+1359>:	movaps 0x40(%rcx),%xmm13
   0x0000000000000563 <+1379>:	movaps (%r8),%xmm4
   0x0000000000000571 <+1393>:	movaps 0x40(%r8),%xmm14
   0x0000000000000585 <+1413>:	mulps  %xmm13,%xmm4
   0x00000000000005ae <+1454>:	movaps 0x80(%r8),%xmm12
   0x00000000000005b6 <+1462>:	mulps  %xmm13,%xmm14
   0x00000000000005ea <+1514>:	mulps  %xmm13,%xmm12
   0x00000000000005f6 <+1526>:	movaps 0xc0(%r8),%xmm12
   0x000000000000061d <+1565>:	mulps  %xmm13,%xmm12

221
222	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
223	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
224	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000554 <+1364>:	movaps 0x50(%rcx),%xmm8
   0x0000000000000567 <+1383>:	movaps 0x10(%r8),%xmm9
   0x0000000000000576 <+1398>:	movaps 0x50(%r8),%xmm0
   0x0000000000000589 <+1417>:	mulps  %xmm8,%xmm9
   0x000000000000058d <+1421>:	addps  %xmm4,%xmm9
   0x00000000000005ba <+1466>:	mulps  %xmm8,%xmm0
   0x00000000000005be <+1470>:	addps  %xmm14,%xmm0
   0x00000000000005c2 <+1474>:	movaps 0x90(%r8),%xmm14
   0x00000000000005ee <+1518>:	mulps  %xmm8,%xmm14
   0x00000000000005f2 <+1522>:	addps  %xmm12,%xmm14
   0x0000000000000621 <+1569>:	movaps 0xd0(%r8),%xmm13
   0x0000000000000629 <+1577>:	mulps  %xmm8,%xmm13
   0x0000000000000635 <+1589>:	addps  %xmm12,%xmm13

225
226	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
227	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
228	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000559 <+1369>:	movaps 0x60(%rcx),%xmm11
   0x000000000000056c <+1388>:	movaps 0x20(%r8),%xmm12
   0x000000000000057b <+1403>:	movaps 0x60(%r8),%xmm15
   0x0000000000000599 <+1433>:	mulps  %xmm11,%xmm12
   0x000000000000059d <+1437>:	addps  %xmm9,%xmm12
   0x00000000000005ca <+1482>:	mulps  %xmm11,%xmm15
   0x00000000000005ce <+1486>:	addps  %xmm0,%xmm15
   0x00000000000005d2 <+1490>:	movaps 0xa0(%r8),%xmm0
   0x00000000000005fe <+1534>:	mulps  %xmm11,%xmm0
   0x0000000000000602 <+1538>:	addps  %xmm14,%xmm0
   0x000000000000062d <+1581>:	movaps 0xe0(%r8),%xmm8
   0x0000000000000639 <+1593>:	mulps  %xmm11,%xmm8
   0x0000000000000645 <+1605>:	addps  %xmm13,%xmm8

229
230	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
231	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
232	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000055e <+1374>:	movaps 0x70(%rcx),%xmm10
   0x0000000000000580 <+1408>:	movaps 0x70(%r8),%xmm3
   0x0000000000000591 <+1425>:	movaps 0xb0(%r8),%xmm4
   0x00000000000005a1 <+1441>:	movaps 0x30(%r8),%xmm9
   0x00000000000005a6 <+1446>:	mulps  %xmm10,%xmm9
   0x00000000000005aa <+1450>:	addps  %xmm12,%xmm9
   0x00000000000005da <+1498>:	mulps  %xmm10,%xmm3
   0x00000000000005de <+1502>:	addps  %xmm15,%xmm3
   0x000000000000060e <+1550>:	mulps  %xmm10,%xmm4
   0x0000000000000612 <+1554>:	addps  %xmm0,%xmm4
   0x0000000000000615 <+1557>:	movaps 0xf0(%r8),%xmm0
   0x0000000000000651 <+1617>:	mulps  %xmm10,%xmm0
   0x0000000000000661 <+1633>:	addps  %xmm8,%xmm0

233	      }
234
235	      /* A(1,2)*B(2,2) = C(1,2). */
236	      for (i = 0; i < 4; i++)
237	      {
238	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
239	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
240	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000606 <+1542>:	movaps 0x140(%rcx),%xmm14
   0x000000000000063d <+1597>:	movaps 0x100(%r8),%xmm11
   0x000000000000065d <+1629>:	mulps  %xmm14,%xmm11
   0x0000000000000671 <+1649>:	addps  %xmm9,%xmm11
   0x0000000000000691 <+1681>:	movaps 0x140(%r8),%xmm11
   0x00000000000006a5 <+1701>:	mulps  %xmm14,%xmm11
   0x00000000000006ad <+1709>:	movaps 0x180(%r8),%xmm13
   0x00000000000006b9 <+1721>:	addps  %xmm3,%xmm11
   0x00000000000006dd <+1757>:	mulps  %xmm14,%xmm13
   0x00000000000006e1 <+1761>:	movaps 0x1c0(%r8),%xmm3
   0x00000000000006f9 <+1785>:	addps  %xmm4,%xmm13
   0x000000000000071d <+1821>:	mulps  %xmm14,%xmm3
   0x0000000000000739 <+1849>:	addps  %xmm0,%xmm3

241
242	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
243	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
244	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000655 <+1621>:	movaps 0x110(%r8),%xmm10
   0x0000000000000675 <+1653>:	movaps 0x150(%rcx),%xmm9
   0x0000000000000685 <+1669>:	mulps  %xmm9,%xmm10
   0x0000000000000689 <+1673>:	addps  %xmm11,%xmm10
   0x00000000000006bd <+1725>:	movaps 0x150(%r8),%xmm3
   0x00000000000006c5 <+1733>:	mulps  %xmm9,%xmm3
   0x00000000000006c9 <+1737>:	addps  %xmm11,%xmm3
   0x00000000000006fd <+1789>:	movaps 0x190(%r8),%xmm4
   0x0000000000000705 <+1797>:	mulps  %xmm9,%xmm4
   0x0000000000000709 <+1801>:	addps  %xmm13,%xmm4
   0x0000000000000721 <+1825>:	movaps 0x1d0(%r8),%xmm14
   0x000000000000072d <+1837>:	mulps  %xmm9,%xmm14
   0x000000000000075c <+1884>:	addps  %xmm3,%xmm14

245
246	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
247	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
248	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000005e2 <+1506>:	movaps 0x160(%rcx),%xmm15
   0x0000000000000649 <+1609>:	movaps 0x120(%r8),%xmm13
   0x000000000000066d <+1645>:	mulps  %xmm15,%xmm13
   0x0000000000000699 <+1689>:	addps  %xmm10,%xmm13
   0x00000000000006cd <+1741>:	movaps 0x160(%r8),%xmm11
   0x00000000000006d5 <+1749>:	mulps  %xmm15,%xmm11
   0x00000000000006d9 <+1753>:	addps  %xmm3,%xmm11
   0x000000000000070d <+1805>:	movaps 0x1a0(%r8),%xmm13
   0x0000000000000715 <+1813>:	mulps  %xmm15,%xmm13
   0x0000000000000719 <+1817>:	addps  %xmm4,%xmm13
   0x000000000000073c <+1852>:	movaps 0x1e0(%r8),%xmm0
   0x0000000000000744 <+1860>:	mulps  %xmm15,%xmm0
   0x000000000000077a <+1914>:	addps  %xmm14,%xmm0

249
250	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
251	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
252	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000665 <+1637>:	movaps 0x130(%r8),%xmm8
   0x000000000000067d <+1661>:	movaps 0x170(%rcx),%xmm12
   0x000000000000068d <+1677>:	mulps  %xmm12,%xmm8
   0x000000000000069d <+1693>:	movaps 0x170(%r8),%xmm10
   0x00000000000006a9 <+1705>:	addps  %xmm13,%xmm8
   0x00000000000006b5 <+1717>:	mulps  %xmm12,%xmm10
   0x00000000000006e9 <+1769>:	addps  %xmm11,%xmm10
   0x00000000000006ed <+1773>:	movaps 0x1b0(%r8),%xmm11
   0x00000000000006f5 <+1781>:	mulps  %xmm12,%xmm11
   0x0000000000000729 <+1833>:	addps  %xmm13,%xmm11
   0x0000000000000731 <+1841>:	movaps 0x1f0(%r8),%xmm9
   0x0000000000000750 <+1872>:	mulps  %xmm12,%xmm9
   0x000000000000078a <+1930>:	addps  %xmm0,%xmm9

253	      }
254
255	      /* A(1,3)*B(3,2) = C(1,2). */
256	      for (i = 0; i < 4; i++)
257	      {
258	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
259	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
260	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000754 <+1876>:	movaps 0x200(%r8),%xmm12
   0x000000000000078e <+1934>:	movaps 0x240(%rcx),%xmm0
   0x0000000000000795 <+1941>:	mulps  %xmm0,%xmm12
   0x0000000000000799 <+1945>:	addps  %xmm8,%xmm12
   0x000000000000079d <+1949>:	movaps 0x240(%r8),%xmm8
   0x00000000000007a5 <+1957>:	mulps  %xmm0,%xmm8
   0x00000000000007bd <+1981>:	movaps 0x280(%r8),%xmm15
   0x00000000000007c5 <+1989>:	mulps  %xmm0,%xmm15
   0x00000000000007d9 <+2009>:	addps  %xmm10,%xmm8
   0x00000000000007ed <+2029>:	movaps 0x2c0(%r8),%xmm8
   0x00000000000007f5 <+2037>:	mulps  %xmm0,%xmm8
   0x0000000000000809 <+2057>:	addps  %xmm11,%xmm15
   0x0000000000000869 <+2153>:	addps  %xmm9,%xmm8

261
262	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
263	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
264	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000748 <+1864>:	movaps 0x210(%r8),%xmm15
   0x0000000000000760 <+1888>:	movaps 0x250(%rcx),%xmm3
   0x0000000000000776 <+1910>:	mulps  %xmm3,%xmm15
   0x00000000000007a9 <+1961>:	addps  %xmm12,%xmm15
   0x00000000000007cd <+1997>:	movaps 0x250(%r8),%xmm14
   0x00000000000007d5 <+2005>:	mulps  %xmm3,%xmm14
   0x00000000000007e9 <+2025>:	addps  %xmm8,%xmm14
   0x000000000000080d <+2061>:	movaps 0x290(%r8),%xmm11
   0x0000000000000815 <+2069>:	mulps  %xmm3,%xmm11
   0x0000000000000829 <+2089>:	addps  %xmm15,%xmm11
   0x000000000000086d <+2157>:	movaps 0x2d0(%r8),%xmm9
   0x0000000000000875 <+2165>:	mulps  %xmm3,%xmm9
   0x0000000000000879 <+2169>:	addps  %xmm8,%xmm9

265
266	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
267	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
268	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000767 <+1895>:	movaps 0x260(%rcx),%xmm13
   0x000000000000077e <+1918>:	movaps 0x220(%r8),%xmm14
   0x0000000000000786 <+1926>:	mulps  %xmm13,%xmm14
   0x00000000000007b9 <+1977>:	addps  %xmm15,%xmm14
   0x00000000000007dd <+2013>:	movaps 0x260(%r8),%xmm10
   0x00000000000007e5 <+2021>:	mulps  %xmm13,%xmm10
   0x00000000000007f9 <+2041>:	addps  %xmm14,%xmm10
   0x000000000000081d <+2077>:	movaps 0x2a0(%r8),%xmm10
   0x0000000000000825 <+2085>:	mulps  %xmm13,%xmm10
   0x0000000000000839 <+2105>:	addps  %xmm11,%xmm10
   0x000000000000083d <+2109>:	movaps 0x2e0(%r8),%xmm11
   0x0000000000000845 <+2117>:	mulps  %xmm13,%xmm11
   0x0000000000000890 <+2192>:	addps  %xmm9,%xmm11

269
270	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
271	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
272	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000076f <+1903>:	movaps 0x270(%rcx),%xmm4
   0x00000000000007ad <+1965>:	movaps 0x230(%r8),%xmm12
   0x00000000000007b5 <+1973>:	mulps  %xmm4,%xmm12
   0x00000000000007c9 <+1993>:	addps  %xmm14,%xmm12
   0x00000000000007fd <+2045>:	movaps 0x270(%r8),%xmm14
   0x0000000000000805 <+2053>:	mulps  %xmm4,%xmm14
   0x0000000000000819 <+2073>:	addps  %xmm10,%xmm14
   0x000000000000082d <+2093>:	movaps 0x2b0(%r8),%xmm15
   0x0000000000000835 <+2101>:	mulps  %xmm4,%xmm15
   0x0000000000000849 <+2121>:	movaps 0x2f0(%r8),%xmm13
   0x0000000000000851 <+2129>:	mulps  %xmm4,%xmm13
   0x0000000000000855 <+2133>:	addps  %xmm10,%xmm15
   0x000000000000089f <+2207>:	addps  %xmm11,%xmm13

273	      }
274
275	      /* A(1,4)*B(4,2) = C(1,2). */
276	      for (i = 0; i < 4; i++)
277	      {
278	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
279	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
280	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000859 <+2137>:	movaps 0x300(%r8),%xmm10
   0x0000000000000861 <+2145>:	movaps 0x340(%r8),%xmm4
   0x000000000000087d <+2173>:	movaps 0x340(%rcx),%xmm0
   0x000000000000088c <+2188>:	mulps  %xmm0,%xmm10
   0x000000000000089c <+2204>:	mulps  %xmm0,%xmm4
   0x00000000000008b7 <+2231>:	addps  %xmm12,%xmm10
   0x00000000000008c7 <+2247>:	addps  %xmm14,%xmm4
   0x00000000000008f7 <+2295>:	movaps 0x380(%r8),%xmm3
   0x0000000000000907 <+2311>:	movaps 0x3c0(%r8),%xmm4
   0x000000000000090f <+2319>:	mulps  %xmm0,%xmm3
   0x0000000000000916 <+2326>:	mulps  %xmm0,%xmm4
   0x0000000000000921 <+2337>:	addps  %xmm15,%xmm3
   0x0000000000000939 <+2361>:	addps  %xmm13,%xmm4

281
282	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
283	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
284	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000884 <+2180>:	movaps 0x350(%rcx),%xmm8
   0x00000000000008bb <+2235>:	movaps 0x310(%r8),%xmm12
   0x00000000000008c3 <+2243>:	mulps  %xmm8,%xmm12
   0x00000000000008d3 <+2259>:	addps  %xmm10,%xmm12
   0x00000000000008eb <+2283>:	movaps 0x350(%r8),%xmm12
   0x00000000000008ff <+2303>:	mulps  %xmm8,%xmm12
   0x0000000000000903 <+2307>:	addps  %xmm4,%xmm12
   0x0000000000000919 <+2329>:	movaps 0x3d0(%r8),%xmm0
   0x0000000000000925 <+2341>:	movaps 0x390(%r8),%xmm15
   0x000000000000093d <+2365>:	mulps  %xmm8,%xmm15
   0x0000000000000959 <+2393>:	addps  %xmm3,%xmm15
   0x000000000000096d <+2413>:	mulps  %xmm8,%xmm0
   0x0000000000000979 <+2425>:	addps  %xmm4,%xmm0

285
286	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
287	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
288	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000008a3 <+2211>:	movaps 0x320(%r8),%xmm3
   0x00000000000008cb <+2251>:	movaps 0x360(%r8),%xmm14
   0x00000000000008d7 <+2263>:	movaps 0x360(%rcx),%xmm10
   0x00000000000008df <+2271>:	mulps  %xmm10,%xmm3
   0x00000000000008e3 <+2275>:	addps  %xmm12,%xmm3
   0x00000000000008e7 <+2279>:	mulps  %xmm10,%xmm14
   0x0000000000000912 <+2322>:	addps  %xmm12,%xmm14
   0x0000000000000941 <+2369>:	movaps 0x3e0(%r8),%xmm13
   0x000000000000095d <+2397>:	movaps 0x3a0(%r8),%xmm3
   0x0000000000000965 <+2405>:	mulps  %xmm10,%xmm3
   0x0000000000000969 <+2409>:	addps  %xmm15,%xmm3
   0x0000000000000975 <+2421>:	mulps  %xmm10,%xmm13
   0x00000000000009a6 <+2470>:	addps  %xmm0,%xmm13

289
290	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
291	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
292	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000894 <+2196>:	movaps 0x370(%rcx),%xmm9
   0x00000000000008ab <+2219>:	movaps 0x330(%r8),%xmm11
   0x00000000000008b3 <+2227>:	mulps  %xmm9,%xmm11
   0x00000000000008f3 <+2291>:	addps  %xmm3,%xmm11
   0x000000000000092d <+2349>:	movaps 0x370(%r8),%xmm12
   0x0000000000000935 <+2357>:	mulps  %xmm9,%xmm12
   0x0000000000000949 <+2377>:	addps  %xmm14,%xmm12
   0x000000000000094d <+2381>:	movaps 0x3b0(%r8),%xmm14
   0x0000000000000955 <+2389>:	mulps  %xmm9,%xmm14
   0x0000000000000971 <+2417>:	addps  %xmm3,%xmm14
   0x00000000000009aa <+2474>:	movaps 0x3f0(%r8),%xmm0
   0x00000000000009b2 <+2482>:	mulps  %xmm9,%xmm0
   0x00000000000009b6 <+2486>:	addps  %xmm13,%xmm0

293	      }
294
295	      /* Store C(1,2) block. */
296	      for (i = 0; i < 4; i++)
297	      {
298	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x000000000000097c <+2428>:	mulps  %xmm2,%xmm11
   0x000000000000098a <+2442>:	mulps  %xmm2,%xmm12
   0x0000000000000998 <+2456>:	mulps  %xmm2,%xmm14
   0x00000000000009ba <+2490>:	mulps  %xmm2,%xmm0

299	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_12]), C_row[i]);
   0x0000000000000980 <+2432>:	addps  0x40(%r9),%xmm11
   0x000000000000098e <+2446>:	addps  0x50(%r9),%xmm12
   0x000000000000099c <+2460>:	addps  0x60(%r9),%xmm14
   0x00000000000009bd <+2493>:	addps  0x70(%r9),%xmm0

300	        _mm_store_ps(&C[i*4+C_OFFSET_12], C_row[i]);
   0x0000000000000985 <+2437>:	movaps %xmm11,0x40(%r9)
   0x0000000000000993 <+2451>:	movaps %xmm12,0x50(%r9)
   0x00000000000009a1 <+2465>:	movaps %xmm14,0x60(%r9)
   0x00000000000009c2 <+2498>:	movaps %xmm0,0x70(%r9)

301	      }
302	    }
303
304	    /* Reset C(1,3) matrix accumulators */
305	    C_row[0] = _mm_setzero_ps();
306	    C_row[1] = _mm_setzero_ps();
307	    C_row[2] = _mm_setzero_ps();
308	    C_row[3] = _mm_setzero_ps();
309
310	    if (norm[0]*norm[18] >= tolerance &&
   0x00000000000009c7 <+2503>:	movss  0x60(%rdx,%rsi,1),%xmm4
   0x00000000000009cd <+2509>:	movaps %xmm7,%xmm0
   0x00000000000009d0 <+2512>:	mulss  %xmm4,%xmm0
   0x00000000000009d4 <+2516>:	comiss %xmm1,%xmm0
   0x00000000000009d7 <+2519>:	jb     0xec0 <stream_kernel+3776>

311	        norm[1]*norm[22] >= tolerance &&
   0x00000000000009dd <+2525>:	movss  0x1c(%rdx,%rsi,1),%xmm0
   0x00000000000009e3 <+2531>:	mulss  0x70(%rdx,%rsi,1),%xmm0
   0x00000000000009e9 <+2537>:	comiss %xmm1,%xmm0
   0x00000000000009ec <+2540>:	jb     0xec0 <stream_kernel+3776>

312	        norm[2]*norm[26] >= tolerance &&
   0x00000000000009f2 <+2546>:	movss  0x20(%rdx,%rsi,1),%xmm0
   0x00000000000009f8 <+2552>:	mulss  0x80(%rdx,%rsi,1),%xmm0
   0x0000000000000a01 <+2561>:	comiss %xmm1,%xmm0
   0x0000000000000a04 <+2564>:	jb     0xec0 <stream_kernel+3776>

313	        norm[3]*norm[30] >= tolerance)
   0x0000000000000a0a <+2570>:	movss  0x24(%rdx,%rsi,1),%xmm0
   0x0000000000000a10 <+2576>:	mulss  0x90(%rdx,%rsi,1),%xmm0
   0x0000000000000a19 <+2585>:	comiss %xmm1,%xmm0
   0x0000000000000a1c <+2588>:	jb     0xec0 <stream_kernel+3776>

314	    {
315	      /* A(1,1)*B(1,3) = C(1,3). */
316	      for (i = 0; i < 4; i++)
317	      {
318	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
319	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
320	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000a22 <+2594>:	movaps 0x80(%rcx),%xmm3
   0x0000000000000a40 <+2624>:	movaps (%r8),%xmm10
   0x0000000000000a53 <+2643>:	movaps 0x40(%r8),%xmm15
   0x0000000000000a5d <+2653>:	mulps  %xmm3,%xmm10
   0x0000000000000a83 <+2691>:	movaps 0x80(%r8),%xmm11
   0x0000000000000a8b <+2699>:	mulps  %xmm3,%xmm15
   0x0000000000000aa7 <+2727>:	movaps 0xc0(%r8),%xmm14
   0x0000000000000abf <+2751>:	mulps  %xmm3,%xmm11
   0x0000000000000af3 <+2803>:	mulps  %xmm3,%xmm14

321
322	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
323	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
324	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000a29 <+2601>:	movaps 0x90(%rcx),%xmm0
   0x0000000000000a44 <+2628>:	movaps 0x10(%r8),%xmm9
   0x0000000000000a58 <+2648>:	movaps 0x50(%r8),%xmm14
   0x0000000000000a61 <+2657>:	mulps  %xmm0,%xmm9
   0x0000000000000a65 <+2661>:	addps  %xmm10,%xmm9
   0x0000000000000a8f <+2703>:	mulps  %xmm0,%xmm14
   0x0000000000000a93 <+2707>:	addps  %xmm15,%xmm14
   0x0000000000000a97 <+2711>:	movaps 0x90(%r8),%xmm15
   0x0000000000000ac3 <+2755>:	mulps  %xmm0,%xmm15
   0x0000000000000ac7 <+2759>:	addps  %xmm11,%xmm15
   0x0000000000000af7 <+2807>:	movaps 0xd0(%r8),%xmm3
   0x0000000000000aff <+2815>:	mulps  %xmm0,%xmm3
   0x0000000000000b0a <+2826>:	addps  %xmm14,%xmm3

325
326	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
327	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
328	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000a30 <+2608>:	movaps 0xa0(%rcx),%xmm8
   0x0000000000000a49 <+2633>:	movaps 0x20(%r8),%xmm11
   0x0000000000000a69 <+2665>:	movaps 0x60(%r8),%xmm10
   0x0000000000000a6e <+2670>:	mulps  %xmm8,%xmm11
   0x0000000000000a72 <+2674>:	addps  %xmm9,%xmm11
   0x0000000000000a9f <+2719>:	mulps  %xmm8,%xmm10
   0x0000000000000aa3 <+2723>:	addps  %xmm14,%xmm10
   0x0000000000000acb <+2763>:	movaps 0xa0(%r8),%xmm11
   0x0000000000000ad3 <+2771>:	mulps  %xmm8,%xmm11
   0x0000000000000ad7 <+2775>:	addps  %xmm15,%xmm11
   0x0000000000000b02 <+2818>:	movaps 0xe0(%r8),%xmm0
   0x0000000000000b0e <+2830>:	mulps  %xmm8,%xmm0
   0x0000000000000b26 <+2854>:	addps  %xmm3,%xmm0

329
330	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
331	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
332	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000a38 <+2616>:	movaps 0xb0(%rcx),%xmm12
   0x0000000000000a4e <+2638>:	movaps 0x30(%r8),%xmm13
   0x0000000000000a76 <+2678>:	movaps 0x70(%r8),%xmm9
   0x0000000000000a7b <+2683>:	mulps  %xmm12,%xmm13
   0x0000000000000a7f <+2687>:	addps  %xmm11,%xmm13
   0x0000000000000aaf <+2735>:	mulps  %xmm12,%xmm9
   0x0000000000000ab3 <+2739>:	addps  %xmm10,%xmm9
   0x0000000000000ab7 <+2743>:	movaps 0xb0(%r8),%xmm10
   0x0000000000000ae3 <+2787>:	mulps  %xmm12,%xmm10
   0x0000000000000ae7 <+2791>:	addps  %xmm11,%xmm10
   0x0000000000000aeb <+2795>:	movaps 0xf0(%r8),%xmm11
   0x0000000000000b22 <+2850>:	mulps  %xmm12,%xmm11
   0x0000000000000b3d <+2877>:	addps  %xmm0,%xmm11

333	      }
334
335	      /* A(1,2)*B(2,3) = C(1,3). */
336	      for (i = 0; i < 4; i++)
337	      {
338	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
339	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
340	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000b12 <+2834>:	movaps 0x180(%rcx),%xmm8
   0x0000000000000b29 <+2857>:	movaps 0x100(%r8),%xmm12
   0x0000000000000b39 <+2873>:	mulps  %xmm8,%xmm12
   0x0000000000000b41 <+2881>:	movaps 0x140(%r8),%xmm0
   0x0000000000000b4d <+2893>:	addps  %xmm13,%xmm12
   0x0000000000000b61 <+2913>:	mulps  %xmm8,%xmm0
   0x0000000000000b7d <+2941>:	addps  %xmm9,%xmm0
   0x0000000000000bb1 <+2993>:	movaps 0x180(%r8),%xmm9
   0x0000000000000bb9 <+3001>:	mulps  %xmm8,%xmm9
   0x0000000000000bcc <+3020>:	addps  %xmm10,%xmm9
   0x0000000000000bf0 <+3056>:	movaps 0x1c0(%r8),%xmm10
   0x0000000000000bf8 <+3064>:	mulps  %xmm8,%xmm10
   0x0000000000000c1c <+3100>:	addps  %xmm11,%xmm10

341
342	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
343	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
344	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000b1a <+2842>:	movaps 0x190(%rcx),%xmm14
   0x0000000000000b51 <+2897>:	movaps 0x110(%r8),%xmm13
   0x0000000000000b59 <+2905>:	mulps  %xmm14,%xmm13
   0x0000000000000b5d <+2909>:	addps  %xmm12,%xmm13
   0x0000000000000b81 <+2945>:	movaps 0x150(%r8),%xmm9
   0x0000000000000b89 <+2953>:	mulps  %xmm14,%xmm9
   0x0000000000000b9d <+2973>:	addps  %xmm0,%xmm9
   0x0000000000000bd0 <+3024>:	movaps 0x190(%r8),%xmm10
   0x0000000000000bd8 <+3032>:	mulps  %xmm14,%xmm10
   0x0000000000000bdc <+3036>:	addps  %xmm9,%xmm10
   0x0000000000000bfc <+3068>:	movaps 0x1d0(%r8),%xmm8
   0x0000000000000c04 <+3076>:	mulps  %xmm14,%xmm8
   0x0000000000000c38 <+3128>:	addps  %xmm10,%xmm8

345
346	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
347	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
348	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000adb <+2779>:	movaps 0x1a0(%rcx),%xmm15
   0x0000000000000b31 <+2865>:	movaps 0x120(%r8),%xmm3
   0x0000000000000b49 <+2889>:	mulps  %xmm15,%xmm3
   0x0000000000000b6d <+2925>:	addps  %xmm13,%xmm3
   0x0000000000000ba1 <+2977>:	movaps 0x160(%r8),%xmm0
   0x0000000000000ba9 <+2985>:	mulps  %xmm15,%xmm0
   0x0000000000000bad <+2989>:	addps  %xmm9,%xmm0
   0x0000000000000be0 <+3040>:	movaps 0x1a0(%r8),%xmm9
   0x0000000000000be8 <+3048>:	mulps  %xmm15,%xmm9
   0x0000000000000bec <+3052>:	addps  %xmm10,%xmm9
   0x0000000000000c3c <+3132>:	movaps 0x1e0(%r8),%xmm10
   0x0000000000000c44 <+3140>:	mulps  %xmm15,%xmm10
   0x0000000000000c58 <+3160>:	addps  %xmm8,%xmm10

349
350	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
351	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
352	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000b65 <+2917>:	movaps 0x1b0(%rcx),%xmm12
   0x0000000000000b71 <+2929>:	movaps 0x130(%r8),%xmm13
   0x0000000000000b79 <+2937>:	mulps  %xmm12,%xmm13
   0x0000000000000b8d <+2957>:	addps  %xmm3,%xmm13
   0x0000000000000b91 <+2961>:	movaps 0x170(%r8),%xmm3
   0x0000000000000b99 <+2969>:	mulps  %xmm12,%xmm3
   0x0000000000000bbd <+3005>:	addps  %xmm0,%xmm3
   0x0000000000000bc0 <+3008>:	movaps 0x1b0(%r8),%xmm0
   0x0000000000000bc8 <+3016>:	mulps  %xmm12,%xmm0
   0x0000000000000c08 <+3080>:	addps  %xmm9,%xmm0
   0x0000000000000c20 <+3104>:	movaps 0x1f0(%r8),%xmm11
   0x0000000000000c28 <+3112>:	mulps  %xmm12,%xmm11
   0x0000000000000c64 <+3172>:	addps  %xmm10,%xmm11

353	      }
354
355	      /* A(1,3)*B(3,3) = C(1,3). */
356	      for (i = 0; i < 4; i++)
357	      {
358	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
359	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
360	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000c0c <+3084>:	movaps 0x280(%rcx),%xmm14
   0x0000000000000c2c <+3116>:	movaps 0x200(%r8),%xmm12
   0x0000000000000c34 <+3124>:	mulps  %xmm14,%xmm12
   0x0000000000000c48 <+3144>:	addps  %xmm13,%xmm12
   0x0000000000000c80 <+3200>:	movaps 0x240(%r8),%xmm12
   0x0000000000000c88 <+3208>:	mulps  %xmm14,%xmm12
   0x0000000000000c9c <+3228>:	addps  %xmm3,%xmm12
   0x0000000000000cc0 <+3264>:	movaps 0x280(%r8),%xmm12
   0x0000000000000cc8 <+3272>:	mulps  %xmm14,%xmm12
   0x0000000000000cdc <+3292>:	addps  %xmm0,%xmm12
   0x0000000000000d00 <+3328>:	movaps 0x2c0(%r8),%xmm12
   0x0000000000000d08 <+3336>:	mulps  %xmm14,%xmm12
   0x0000000000000d28 <+3368>:	addps  %xmm11,%xmm12

361
362	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
363	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
364	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000c14 <+3092>:	movaps 0x290(%rcx),%xmm9
   0x0000000000000c4c <+3148>:	movaps 0x210(%r8),%xmm13
   0x0000000000000c54 <+3156>:	mulps  %xmm9,%xmm13
   0x0000000000000c7c <+3196>:	addps  %xmm12,%xmm13
   0x0000000000000ca0 <+3232>:	movaps 0x250(%r8),%xmm3
   0x0000000000000ca8 <+3240>:	mulps  %xmm9,%xmm3
   0x0000000000000cbc <+3260>:	addps  %xmm12,%xmm3
   0x0000000000000ce0 <+3296>:	movaps 0x290(%r8),%xmm0
   0x0000000000000ce8 <+3304>:	mulps  %xmm9,%xmm0
   0x0000000000000cfc <+3324>:	addps  %xmm12,%xmm0
   0x0000000000000d2c <+3372>:	movaps 0x2d0(%r8),%xmm11
   0x0000000000000d34 <+3380>:	mulps  %xmm9,%xmm11
   0x0000000000000d58 <+3416>:	addps  %xmm12,%xmm11

365
366	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
367	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
368	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000c5c <+3164>:	movaps 0x220(%r8),%xmm8
   0x0000000000000c68 <+3176>:	movaps 0x2a0(%rcx),%xmm10
   0x0000000000000c78 <+3192>:	mulps  %xmm10,%xmm8
   0x0000000000000c8c <+3212>:	addps  %xmm13,%xmm8
   0x0000000000000cb0 <+3248>:	movaps 0x260(%r8),%xmm8
   0x0000000000000cb8 <+3256>:	mulps  %xmm10,%xmm8
   0x0000000000000ccc <+3276>:	addps  %xmm3,%xmm8
   0x0000000000000cf0 <+3312>:	movaps 0x2a0(%r8),%xmm8
   0x0000000000000cf8 <+3320>:	mulps  %xmm10,%xmm8
   0x0000000000000d18 <+3352>:	addps  %xmm0,%xmm8
   0x0000000000000d38 <+3384>:	movaps 0x2e0(%r8),%xmm9
   0x0000000000000d40 <+3392>:	mulps  %xmm10,%xmm9
   0x0000000000000d64 <+3428>:	addps  %xmm11,%xmm9

369
370	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
371	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
372	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000c70 <+3184>:	movaps 0x2b0(%rcx),%xmm15
   0x0000000000000c90 <+3216>:	movaps 0x230(%r8),%xmm13
   0x0000000000000c98 <+3224>:	mulps  %xmm15,%xmm13
   0x0000000000000cac <+3244>:	addps  %xmm8,%xmm13
   0x0000000000000cd0 <+3280>:	movaps 0x270(%r8),%xmm3
   0x0000000000000cd8 <+3288>:	mulps  %xmm15,%xmm3
   0x0000000000000cec <+3308>:	addps  %xmm8,%xmm3
   0x0000000000000d0c <+3340>:	movaps 0x2f0(%r8),%xmm14
   0x0000000000000d14 <+3348>:	mulps  %xmm15,%xmm14
   0x0000000000000d1c <+3356>:	movaps 0x2b0(%r8),%xmm0
   0x0000000000000d24 <+3364>:	mulps  %xmm15,%xmm0
   0x0000000000000d4c <+3404>:	addps  %xmm8,%xmm0
   0x0000000000000d7c <+3452>:	addps  %xmm9,%xmm14

373	      }
374
375	      /* A(1,4)*B(4,3) = C(1,3). */
376	      for (i = 0; i < 4; i++)
377	      {
378	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
379	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
380	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d44 <+3396>:	movaps 0x300(%r8),%xmm10
   0x0000000000000d68 <+3432>:	movaps 0x380(%rcx),%xmm15
   0x0000000000000d78 <+3448>:	mulps  %xmm15,%xmm10
   0x0000000000000d8c <+3468>:	addps  %xmm13,%xmm10
   0x0000000000000db0 <+3504>:	movaps 0x340(%r8),%xmm8
   0x0000000000000db8 <+3512>:	mulps  %xmm15,%xmm8
   0x0000000000000dcc <+3532>:	addps  %xmm3,%xmm8
   0x0000000000000df0 <+3568>:	movaps 0x380(%r8),%xmm3
   0x0000000000000df8 <+3576>:	mulps  %xmm15,%xmm3
   0x0000000000000e0c <+3596>:	addps  %xmm0,%xmm3
   0x0000000000000e2d <+3629>:	movaps 0x3c0(%r8),%xmm0
   0x0000000000000e35 <+3637>:	mulps  %xmm15,%xmm0
   0x0000000000000e85 <+3717>:	addps  %xmm14,%xmm0

381
382	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
383	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
384	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d50 <+3408>:	movaps 0x310(%r8),%xmm8
   0x0000000000000d70 <+3440>:	movaps 0x390(%rcx),%xmm11
   0x0000000000000d80 <+3456>:	mulps  %xmm11,%xmm8
   0x0000000000000d9c <+3484>:	addps  %xmm10,%xmm8
   0x0000000000000dd0 <+3536>:	movaps 0x350(%r8),%xmm3
   0x0000000000000dd8 <+3544>:	mulps  %xmm11,%xmm3
   0x0000000000000ddc <+3548>:	addps  %xmm8,%xmm3
   0x0000000000000e0f <+3599>:	movaps 0x390(%r8),%xmm0
   0x0000000000000e17 <+3607>:	mulps  %xmm11,%xmm0
   0x0000000000000e1b <+3611>:	addps  %xmm3,%xmm0
   0x0000000000000e89 <+3721>:	movaps 0x3d0(%r8),%xmm14
   0x0000000000000e91 <+3729>:	mulps  %xmm11,%xmm14
   0x0000000000000e95 <+3733>:	addps  %xmm0,%xmm14

385
386	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
387	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
388	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d5c <+3420>:	movaps 0x320(%r8),%xmm12
   0x0000000000000d90 <+3472>:	movaps 0x3a0(%rcx),%xmm13
   0x0000000000000d98 <+3480>:	mulps  %xmm13,%xmm12
   0x0000000000000dac <+3500>:	addps  %xmm8,%xmm12
   0x0000000000000dc0 <+3520>:	movaps 0x360(%r8),%xmm12
   0x0000000000000dc8 <+3528>:	mulps  %xmm13,%xmm12
   0x0000000000000dec <+3564>:	addps  %xmm3,%xmm12
   0x0000000000000e1e <+3614>:	movaps 0x3a0(%r8),%xmm3
   0x0000000000000e26 <+3622>:	mulps  %xmm13,%xmm3
   0x0000000000000e2a <+3626>:	addps  %xmm0,%xmm3
   0x0000000000000e39 <+3641>:	movaps 0x3e0(%r8),%xmm15
   0x0000000000000e41 <+3649>:	mulps  %xmm13,%xmm15
   0x0000000000000ea5 <+3749>:	addps  %xmm14,%xmm15

389
390	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
391	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
392	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d84 <+3460>:	movaps 0x330(%r8),%xmm9
   0x0000000000000da0 <+3488>:	movaps 0x3b0(%rcx),%xmm10
   0x0000000000000da8 <+3496>:	mulps  %xmm10,%xmm9
   0x0000000000000dbc <+3516>:	addps  %xmm12,%xmm9
   0x0000000000000de0 <+3552>:	movaps 0x370(%r8),%xmm8
   0x0000000000000de8 <+3560>:	mulps  %xmm10,%xmm8
   0x0000000000000dfc <+3580>:	addps  %xmm12,%xmm8
   0x0000000000000e00 <+3584>:	movaps 0x3b0(%r8),%xmm12
   0x0000000000000e08 <+3592>:	mulps  %xmm10,%xmm12
   0x0000000000000e45 <+3653>:	addps  %xmm3,%xmm12
   0x0000000000000e99 <+3737>:	movaps 0x3f0(%r8),%xmm0
   0x0000000000000ea1 <+3745>:	mulps  %xmm10,%xmm0
   0x0000000000000ea9 <+3753>:	addps  %xmm15,%xmm0

393	      }
394
395	      /* Store C(1,3) block. */
396	      for (i = 0; i < 4; i++)
397	      {
398	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000000e49 <+3657>:	mulps  %xmm2,%xmm9
   0x0000000000000e5d <+3677>:	mulps  %xmm2,%xmm8
   0x0000000000000e71 <+3697>:	mulps  %xmm2,%xmm12
   0x0000000000000ead <+3757>:	mulps  %xmm2,%xmm0

399	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_13]), C_row[i]);
   0x0000000000000e4d <+3661>:	addps  0x80(%r9),%xmm9
   0x0000000000000e61 <+3681>:	addps  0x90(%r9),%xmm8
   0x0000000000000e75 <+3701>:	addps  0xa0(%r9),%xmm12
   0x0000000000000eb0 <+3760>:	addps  0xb0(%r9),%xmm0

400	        _mm_store_ps(&C[i*4+C_OFFSET_13], C_row[i]);
   0x0000000000000e55 <+3669>:	movaps %xmm9,0x80(%r9)
   0x0000000000000e69 <+3689>:	movaps %xmm8,0x90(%r9)
   0x0000000000000e7d <+3709>:	movaps %xmm12,0xa0(%r9)
   0x0000000000000eb8 <+3768>:	movaps %xmm0,0xb0(%r9)

401	      }
402	    }
403
404	    /* Reset C(1,4) matrix accumulators */
405	    C_row[0] = _mm_setzero_ps();
406	    C_row[1] = _mm_setzero_ps();
407	    C_row[2] = _mm_setzero_ps();
408	    C_row[3] = _mm_setzero_ps();
409
410	    if (norm[0]*norm[19] >= tolerance &&
   0x0000000000000ec0 <+3776>:	movss  0x64(%rdx,%rsi,1),%xmm3
   0x0000000000000ec6 <+3782>:	mulss  %xmm3,%xmm7
   0x0000000000000eca <+3786>:	comiss %xmm1,%xmm7
   0x0000000000000ecd <+3789>:	jb     0x13b6 <stream_kernel+5046>

411	        norm[1]*norm[23] >= tolerance &&
   0x0000000000000ed3 <+3795>:	movss  0x1c(%rdx,%rsi,1),%xmm0
   0x0000000000000ed9 <+3801>:	mulss  0x74(%rdx,%rsi,1),%xmm0
   0x0000000000000edf <+3807>:	comiss %xmm1,%xmm0
   0x0000000000000ee2 <+3810>:	jb     0x13b6 <stream_kernel+5046>

412	        norm[2]*norm[27] >= tolerance &&
   0x0000000000000ee8 <+3816>:	movss  0x20(%rdx,%rsi,1),%xmm0
   0x0000000000000eee <+3822>:	mulss  0x84(%rdx,%rsi,1),%xmm0
   0x0000000000000ef7 <+3831>:	comiss %xmm1,%xmm0
   0x0000000000000efa <+3834>:	jb     0x13b6 <stream_kernel+5046>

413	        norm[3]*norm[31] >= tolerance)
   0x0000000000000f00 <+3840>:	movss  0x24(%rdx,%rsi,1),%xmm0
   0x0000000000000f06 <+3846>:	mulss  0x94(%rdx,%rsi,1),%xmm0
   0x0000000000000f0f <+3855>:	comiss %xmm1,%xmm0
   0x0000000000000f12 <+3858>:	jb     0x13b6 <stream_kernel+5046>

414	    {
415	      /* A(1,1)*B(1,4) = C(1,4). */
416	      for (i = 0; i < 4; i++)
417	      {
418	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
419	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
420	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000f18 <+3864>:	movaps 0xc0(%rcx),%xmm7
   0x0000000000000f36 <+3894>:	movaps (%r8),%xmm10
   0x0000000000000f49 <+3913>:	movaps 0x40(%r8),%xmm15
   0x0000000000000f53 <+3923>:	mulps  %xmm7,%xmm10
   0x0000000000000f79 <+3961>:	movaps 0x80(%r8),%xmm11
   0x0000000000000f81 <+3969>:	mulps  %xmm7,%xmm15
   0x0000000000000f9d <+3997>:	movaps 0xc0(%r8),%xmm14
   0x0000000000000fb5 <+4021>:	mulps  %xmm7,%xmm11
   0x0000000000000fe9 <+4073>:	mulps  %xmm7,%xmm14

421
422	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
423	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
424	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000f1f <+3871>:	movaps 0xd0(%rcx),%xmm0
   0x0000000000000f3a <+3898>:	movaps 0x10(%r8),%xmm9
   0x0000000000000f4e <+3918>:	movaps 0x50(%r8),%xmm14
   0x0000000000000f57 <+3927>:	mulps  %xmm0,%xmm9
   0x0000000000000f5b <+3931>:	addps  %xmm10,%xmm9
   0x0000000000000f85 <+3973>:	mulps  %xmm0,%xmm14
   0x0000000000000f89 <+3977>:	addps  %xmm15,%xmm14
   0x0000000000000f8d <+3981>:	movaps 0x90(%r8),%xmm15
   0x0000000000000fb9 <+4025>:	mulps  %xmm0,%xmm15
   0x0000000000000fbd <+4029>:	addps  %xmm11,%xmm15
   0x0000000000000fed <+4077>:	movaps 0xd0(%r8),%xmm7
   0x0000000000000ff5 <+4085>:	mulps  %xmm0,%xmm7
   0x0000000000001000 <+4096>:	addps  %xmm14,%xmm7

425
426	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
427	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
428	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000f26 <+3878>:	movaps 0xe0(%rcx),%xmm8
   0x0000000000000f3f <+3903>:	movaps 0x20(%r8),%xmm11
   0x0000000000000f5f <+3935>:	movaps 0x60(%r8),%xmm10
   0x0000000000000f64 <+3940>:	mulps  %xmm8,%xmm11
   0x0000000000000f68 <+3944>:	addps  %xmm9,%xmm11
   0x0000000000000f95 <+3989>:	mulps  %xmm8,%xmm10
   0x0000000000000f99 <+3993>:	addps  %xmm14,%xmm10
   0x0000000000000fc1 <+4033>:	movaps 0xa0(%r8),%xmm11
   0x0000000000000fc9 <+4041>:	mulps  %xmm8,%xmm11
   0x0000000000000fcd <+4045>:	addps  %xmm15,%xmm11
   0x0000000000000ff8 <+4088>:	movaps 0xe0(%r8),%xmm0
   0x0000000000001004 <+4100>:	mulps  %xmm8,%xmm0
   0x000000000000101c <+4124>:	addps  %xmm7,%xmm0

429
430	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
431	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
432	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000f2e <+3886>:	movaps 0xf0(%rcx),%xmm12
   0x0000000000000f44 <+3908>:	movaps 0x30(%r8),%xmm13
   0x0000000000000f6c <+3948>:	movaps 0x70(%r8),%xmm9
   0x0000000000000f71 <+3953>:	mulps  %xmm12,%xmm13
   0x0000000000000f75 <+3957>:	addps  %xmm11,%xmm13
   0x0000000000000fa5 <+4005>:	mulps  %xmm12,%xmm9
   0x0000000000000fa9 <+4009>:	addps  %xmm10,%xmm9
   0x0000000000000fad <+4013>:	movaps 0xb0(%r8),%xmm10
   0x0000000000000fd9 <+4057>:	mulps  %xmm12,%xmm10
   0x0000000000000fdd <+4061>:	addps  %xmm11,%xmm10
   0x0000000000000fe1 <+4065>:	movaps 0xf0(%r8),%xmm11
   0x0000000000001018 <+4120>:	mulps  %xmm12,%xmm11
   0x0000000000001033 <+4147>:	addps  %xmm0,%xmm11

433	      }
434
435	      /* A(1,2)*B(2,4) = C(1,4). */
436	      for (i = 0; i < 4; i++)
437	      {
438	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
439	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
440	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001008 <+4104>:	movaps 0x1c0(%rcx),%xmm8
   0x000000000000101f <+4127>:	movaps 0x100(%r8),%xmm12
   0x000000000000102f <+4143>:	mulps  %xmm8,%xmm12
   0x0000000000001037 <+4151>:	movaps 0x140(%r8),%xmm0
   0x0000000000001043 <+4163>:	addps  %xmm13,%xmm12
   0x0000000000001057 <+4183>:	mulps  %xmm8,%xmm0
   0x0000000000001073 <+4211>:	addps  %xmm9,%xmm0
   0x00000000000010a7 <+4263>:	movaps 0x180(%r8),%xmm9
   0x00000000000010af <+4271>:	mulps  %xmm8,%xmm9
   0x00000000000010c2 <+4290>:	addps  %xmm10,%xmm9
   0x00000000000010e6 <+4326>:	movaps 0x1c0(%r8),%xmm10
   0x00000000000010ee <+4334>:	mulps  %xmm8,%xmm10
   0x0000000000001112 <+4370>:	addps  %xmm11,%xmm10

441
442	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
443	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
444	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001010 <+4112>:	movaps 0x1d0(%rcx),%xmm14
   0x0000000000001047 <+4167>:	movaps 0x110(%r8),%xmm13
   0x000000000000104f <+4175>:	mulps  %xmm14,%xmm13
   0x0000000000001053 <+4179>:	addps  %xmm12,%xmm13
   0x0000000000001077 <+4215>:	movaps 0x150(%r8),%xmm9
   0x000000000000107f <+4223>:	mulps  %xmm14,%xmm9
   0x0000000000001093 <+4243>:	addps  %xmm0,%xmm9
   0x00000000000010c6 <+4294>:	movaps 0x190(%r8),%xmm10
   0x00000000000010ce <+4302>:	mulps  %xmm14,%xmm10
   0x00000000000010d2 <+4306>:	addps  %xmm9,%xmm10
   0x00000000000010f2 <+4338>:	movaps 0x1d0(%r8),%xmm8
   0x00000000000010fa <+4346>:	mulps  %xmm14,%xmm8
   0x000000000000112e <+4398>:	addps  %xmm10,%xmm8

445
446	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
447	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
448	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000fd1 <+4049>:	movaps 0x1e0(%rcx),%xmm15
   0x0000000000001027 <+4135>:	movaps 0x120(%r8),%xmm7
   0x000000000000103f <+4159>:	mulps  %xmm15,%xmm7
   0x0000000000001063 <+4195>:	addps  %xmm13,%xmm7
   0x0000000000001097 <+4247>:	movaps 0x160(%r8),%xmm0
   0x000000000000109f <+4255>:	mulps  %xmm15,%xmm0
   0x00000000000010a3 <+4259>:	addps  %xmm9,%xmm0
   0x00000000000010d6 <+4310>:	movaps 0x1a0(%r8),%xmm9
   0x00000000000010de <+4318>:	mulps  %xmm15,%xmm9
   0x00000000000010e2 <+4322>:	addps  %xmm10,%xmm9
   0x0000000000001132 <+4402>:	movaps 0x1e0(%r8),%xmm10
   0x000000000000113a <+4410>:	mulps  %xmm15,%xmm10
   0x000000000000114e <+4430>:	addps  %xmm8,%xmm10

449
450	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
451	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
452	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000105b <+4187>:	movaps 0x1f0(%rcx),%xmm12
   0x0000000000001067 <+4199>:	movaps 0x130(%r8),%xmm13
   0x000000000000106f <+4207>:	mulps  %xmm12,%xmm13
   0x0000000000001083 <+4227>:	addps  %xmm7,%xmm13
   0x0000000000001087 <+4231>:	movaps 0x170(%r8),%xmm7
   0x000000000000108f <+4239>:	mulps  %xmm12,%xmm7
   0x00000000000010b3 <+4275>:	addps  %xmm0,%xmm7
   0x00000000000010b6 <+4278>:	movaps 0x1b0(%r8),%xmm0
   0x00000000000010be <+4286>:	mulps  %xmm12,%xmm0
   0x00000000000010fe <+4350>:	addps  %xmm9,%xmm0
   0x0000000000001116 <+4374>:	movaps 0x1f0(%r8),%xmm11
   0x000000000000111e <+4382>:	mulps  %xmm12,%xmm11
   0x000000000000115a <+4442>:	addps  %xmm10,%xmm11

453	      }
454
455	      /* A(1,3)*B(3,4) = C(1,4). */
456	      for (i = 0; i < 4; i++)
457	      {
458	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
459	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
460	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001102 <+4354>:	movaps 0x2c0(%rcx),%xmm14
   0x0000000000001122 <+4386>:	movaps 0x200(%r8),%xmm12
   0x000000000000112a <+4394>:	mulps  %xmm14,%xmm12
   0x000000000000113e <+4414>:	addps  %xmm13,%xmm12
   0x0000000000001176 <+4470>:	movaps 0x240(%r8),%xmm12
   0x000000000000117e <+4478>:	mulps  %xmm14,%xmm12
   0x0000000000001192 <+4498>:	addps  %xmm7,%xmm12
   0x00000000000011b6 <+4534>:	movaps 0x280(%r8),%xmm12
   0x00000000000011be <+4542>:	mulps  %xmm14,%xmm12
   0x00000000000011d2 <+4562>:	addps  %xmm0,%xmm12
   0x00000000000011f6 <+4598>:	movaps 0x2c0(%r8),%xmm12
   0x00000000000011fe <+4606>:	mulps  %xmm14,%xmm12
   0x000000000000121e <+4638>:	addps  %xmm11,%xmm12

461
462	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
463	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
464	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000110a <+4362>:	movaps 0x2d0(%rcx),%xmm9
   0x0000000000001142 <+4418>:	movaps 0x210(%r8),%xmm13
   0x000000000000114a <+4426>:	mulps  %xmm9,%xmm13
   0x0000000000001172 <+4466>:	addps  %xmm12,%xmm13
   0x0000000000001196 <+4502>:	movaps 0x250(%r8),%xmm7
   0x000000000000119e <+4510>:	mulps  %xmm9,%xmm7
   0x00000000000011b2 <+4530>:	addps  %xmm12,%xmm7
   0x00000000000011d6 <+4566>:	movaps 0x290(%r8),%xmm0
   0x00000000000011de <+4574>:	mulps  %xmm9,%xmm0
   0x00000000000011f2 <+4594>:	addps  %xmm12,%xmm0
   0x0000000000001222 <+4642>:	movaps 0x2d0(%r8),%xmm11
   0x000000000000122a <+4650>:	mulps  %xmm9,%xmm11
   0x000000000000124e <+4686>:	addps  %xmm12,%xmm11

465
466	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
467	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
468	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001152 <+4434>:	movaps 0x220(%r8),%xmm8
   0x000000000000115e <+4446>:	movaps 0x2e0(%rcx),%xmm10
   0x000000000000116e <+4462>:	mulps  %xmm10,%xmm8
   0x0000000000001182 <+4482>:	addps  %xmm13,%xmm8
   0x00000000000011a6 <+4518>:	movaps 0x260(%r8),%xmm8
   0x00000000000011ae <+4526>:	mulps  %xmm10,%xmm8
   0x00000000000011c2 <+4546>:	addps  %xmm7,%xmm8
   0x00000000000011e6 <+4582>:	movaps 0x2a0(%r8),%xmm8
   0x00000000000011ee <+4590>:	mulps  %xmm10,%xmm8
   0x000000000000120e <+4622>:	addps  %xmm0,%xmm8
   0x000000000000122e <+4654>:	movaps 0x2e0(%r8),%xmm9
   0x0000000000001236 <+4662>:	mulps  %xmm10,%xmm9
   0x000000000000125a <+4698>:	addps  %xmm11,%xmm9

469
470	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
471	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
472	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001166 <+4454>:	movaps 0x2f0(%rcx),%xmm15
   0x0000000000001186 <+4486>:	movaps 0x230(%r8),%xmm13
   0x000000000000118e <+4494>:	mulps  %xmm15,%xmm13
   0x00000000000011a2 <+4514>:	addps  %xmm8,%xmm13
   0x00000000000011c6 <+4550>:	movaps 0x270(%r8),%xmm7
   0x00000000000011ce <+4558>:	mulps  %xmm15,%xmm7
   0x00000000000011e2 <+4578>:	addps  %xmm8,%xmm7
   0x0000000000001202 <+4610>:	movaps 0x2f0(%r8),%xmm14
   0x000000000000120a <+4618>:	mulps  %xmm15,%xmm14
   0x0000000000001212 <+4626>:	movaps 0x2b0(%r8),%xmm0
   0x000000000000121a <+4634>:	mulps  %xmm15,%xmm0
   0x0000000000001242 <+4674>:	addps  %xmm8,%xmm0
   0x0000000000001272 <+4722>:	addps  %xmm9,%xmm14

473	      }
474
475	      /* A(1,4)*B(4,4) = C(1,4). */
476	      for (i = 0; i < 4; i++)
477	      {
478	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
479	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
480	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000123a <+4666>:	movaps 0x300(%r8),%xmm10
   0x000000000000125e <+4702>:	movaps 0x3c0(%rcx),%xmm15
   0x000000000000126e <+4718>:	mulps  %xmm15,%xmm10
   0x0000000000001282 <+4738>:	addps  %xmm13,%xmm10
   0x00000000000012a6 <+4774>:	movaps 0x340(%r8),%xmm8
   0x00000000000012ae <+4782>:	mulps  %xmm15,%xmm8
   0x00000000000012c2 <+4802>:	addps  %xmm7,%xmm8
   0x00000000000012e6 <+4838>:	movaps 0x380(%r8),%xmm7
   0x00000000000012ee <+4846>:	mulps  %xmm15,%xmm7
   0x0000000000001302 <+4866>:	addps  %xmm0,%xmm7
   0x0000000000001323 <+4899>:	movaps 0x3c0(%r8),%xmm0
   0x000000000000132b <+4907>:	mulps  %xmm15,%xmm0
   0x000000000000137b <+4987>:	addps  %xmm14,%xmm0

481
482	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
483	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
484	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001246 <+4678>:	movaps 0x310(%r8),%xmm8
   0x0000000000001266 <+4710>:	movaps 0x3d0(%rcx),%xmm11
   0x0000000000001276 <+4726>:	mulps  %xmm11,%xmm8
   0x0000000000001292 <+4754>:	addps  %xmm10,%xmm8
   0x00000000000012c6 <+4806>:	movaps 0x350(%r8),%xmm7
   0x00000000000012ce <+4814>:	mulps  %xmm11,%xmm7
   0x00000000000012d2 <+4818>:	addps  %xmm8,%xmm7
   0x0000000000001305 <+4869>:	movaps 0x390(%r8),%xmm0
   0x000000000000130d <+4877>:	mulps  %xmm11,%xmm0
   0x0000000000001311 <+4881>:	addps  %xmm7,%xmm0
   0x000000000000137f <+4991>:	movaps 0x3d0(%r8),%xmm14
   0x0000000000001387 <+4999>:	mulps  %xmm11,%xmm14
   0x000000000000138b <+5003>:	addps  %xmm0,%xmm14

485
486	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
487	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
488	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001252 <+4690>:	movaps 0x320(%r8),%xmm12
   0x0000000000001286 <+4742>:	movaps 0x3e0(%rcx),%xmm13
   0x000000000000128e <+4750>:	mulps  %xmm13,%xmm12
   0x00000000000012a2 <+4770>:	addps  %xmm8,%xmm12
   0x00000000000012b6 <+4790>:	movaps 0x360(%r8),%xmm12
   0x00000000000012be <+4798>:	mulps  %xmm13,%xmm12
   0x00000000000012e2 <+4834>:	addps  %xmm7,%xmm12
   0x0000000000001314 <+4884>:	movaps 0x3a0(%r8),%xmm7
   0x000000000000131c <+4892>:	mulps  %xmm13,%xmm7
   0x0000000000001320 <+4896>:	addps  %xmm0,%xmm7
   0x000000000000132f <+4911>:	movaps 0x3e0(%r8),%xmm15
   0x0000000000001337 <+4919>:	mulps  %xmm13,%xmm15
   0x000000000000139b <+5019>:	addps  %xmm14,%xmm15

489
490	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
491	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
492	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000127a <+4730>:	movaps 0x330(%r8),%xmm9
   0x0000000000001296 <+4758>:	movaps 0x3f0(%rcx),%xmm10
   0x000000000000129e <+4766>:	mulps  %xmm10,%xmm9
   0x00000000000012b2 <+4786>:	addps  %xmm12,%xmm9
   0x00000000000012d6 <+4822>:	movaps 0x370(%r8),%xmm8
   0x00000000000012de <+4830>:	mulps  %xmm10,%xmm8
   0x00000000000012f2 <+4850>:	addps  %xmm12,%xmm8
   0x00000000000012f6 <+4854>:	movaps 0x3b0(%r8),%xmm12
   0x00000000000012fe <+4862>:	mulps  %xmm10,%xmm12
   0x000000000000133b <+4923>:	addps  %xmm7,%xmm12
   0x000000000000138f <+5007>:	movaps 0x3f0(%r8),%xmm0
   0x0000000000001397 <+5015>:	mulps  %xmm10,%xmm0
   0x000000000000139f <+5023>:	addps  %xmm15,%xmm0

493	      }
494
495	      /* Store C(1,4) block. */
496	      for (i = 0; i < 4; i++)
497	      {
498	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x000000000000133f <+4927>:	mulps  %xmm2,%xmm9
   0x0000000000001353 <+4947>:	mulps  %xmm2,%xmm8
   0x0000000000001367 <+4967>:	mulps  %xmm2,%xmm12
   0x00000000000013a3 <+5027>:	mulps  %xmm2,%xmm0

499	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_14]), C_row[i]);
   0x0000000000001343 <+4931>:	addps  0xc0(%r9),%xmm9
   0x0000000000001357 <+4951>:	addps  0xd0(%r9),%xmm8
   0x000000000000136b <+4971>:	addps  0xe0(%r9),%xmm12
   0x00000000000013a6 <+5030>:	addps  0xf0(%r9),%xmm0

500	        _mm_store_ps(&C[i*4+C_OFFSET_14], C_row[i]);
   0x000000000000134b <+4939>:	movaps %xmm9,0xc0(%r9)
   0x000000000000135f <+4959>:	movaps %xmm8,0xd0(%r9)
   0x0000000000001373 <+4979>:	movaps %xmm12,0xe0(%r9)
   0x00000000000013ae <+5038>:	movaps %xmm0,0xf0(%r9)

501	      }
502	    }
503
504	    /* Reset C(2,1) matrix accumulators */
505	    C_row[0] = _mm_setzero_ps();
506	    C_row[1] = _mm_setzero_ps();
507	    C_row[2] = _mm_setzero_ps();
508	    C_row[3] = _mm_setzero_ps();
509
510	    if (norm[4]*norm[16] >= tolerance &&
   0x00000000000013b6 <+5046>:	movss  0x28(%rdx,%rsi,1),%xmm0
   0x00000000000013bc <+5052>:	movaps %xmm6,%xmm7
   0x00000000000013bf <+5055>:	mulss  %xmm0,%xmm7
   0x00000000000013c3 <+5059>:	comiss %xmm1,%xmm7
   0x00000000000013c6 <+5062>:	jb     0x18bf <stream_kernel+6335>

511	        norm[5]*norm[20] >= tolerance &&
   0x00000000000013cc <+5068>:	movss  0x2c(%rdx,%rsi,1),%xmm7
   0x00000000000013d2 <+5074>:	mulss  0x68(%rdx,%rsi,1),%xmm7
   0x00000000000013d8 <+5080>:	comiss %xmm1,%xmm7
   0x00000000000013db <+5083>:	jb     0x18bf <stream_kernel+6335>

512	        norm[6]*norm[24] >= tolerance &&
   0x00000000000013e1 <+5089>:	movss  0x30(%rdx,%rsi,1),%xmm7
   0x00000000000013e7 <+5095>:	mulss  0x78(%rdx,%rsi,1),%xmm7
   0x00000000000013ed <+5101>:	comiss %xmm1,%xmm7
   0x00000000000013f0 <+5104>:	jb     0x18bf <stream_kernel+6335>

513	        norm[7]*norm[28] >= tolerance)
   0x00000000000013f6 <+5110>:	movss  0x34(%rdx,%rsi,1),%xmm7
   0x00000000000013fc <+5116>:	mulss  0x88(%rdx,%rsi,1),%xmm7
   0x0000000000001405 <+5125>:	comiss %xmm1,%xmm7
   0x0000000000001408 <+5128>:	jb     0x18bf <stream_kernel+6335>

514	    {
515	      /* A(2,1)*B(1,1) = C(2,1). */
516	      for (i = 0; i < 4; i++)
517	      {
518	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
519	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
   0x000000000000140e <+5134>:	movaps (%rcx),%xmm7

520	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001420 <+5152>:	movaps 0x400(%r8),%xmm9
   0x0000000000001440 <+5184>:	movaps 0x440(%r8),%xmm15
   0x0000000000001448 <+5192>:	mulps  %xmm7,%xmm9
   0x0000000000001454 <+5204>:	movaps 0x480(%r8),%xmm9
   0x000000000000147c <+5244>:	mulps  %xmm7,%xmm15
   0x00000000000014b0 <+5296>:	mulps  %xmm7,%xmm9
   0x00000000000014bc <+5308>:	movaps 0x4c0(%r8),%xmm9
   0x00000000000014e4 <+5348>:	mulps  %xmm7,%xmm9

521
522	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
523	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
524	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001411 <+5137>:	movaps 0x10(%rcx),%xmm11
   0x0000000000001428 <+5160>:	movaps 0x410(%r8),%xmm13
   0x000000000000144c <+5196>:	mulps  %xmm11,%xmm13
   0x0000000000001450 <+5200>:	addps  %xmm9,%xmm13
   0x0000000000001464 <+5220>:	movaps 0x450(%r8),%xmm13
   0x0000000000001480 <+5248>:	mulps  %xmm11,%xmm13
   0x0000000000001484 <+5252>:	addps  %xmm15,%xmm13
   0x00000000000014a8 <+5288>:	movaps 0x490(%r8),%xmm14
   0x00000000000014b4 <+5300>:	mulps  %xmm11,%xmm14
   0x00000000000014b8 <+5304>:	addps  %xmm9,%xmm14
   0x00000000000014e8 <+5352>:	movaps 0x4d0(%r8),%xmm7
   0x00000000000014f0 <+5360>:	mulps  %xmm11,%xmm7
   0x00000000000014fc <+5372>:	addps  %xmm9,%xmm7

525
526	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
527	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
528	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001416 <+5142>:	movaps 0x20(%rcx),%xmm10
   0x0000000000001430 <+5168>:	movaps 0x420(%r8),%xmm14
   0x000000000000145c <+5212>:	mulps  %xmm10,%xmm14
   0x0000000000001460 <+5216>:	addps  %xmm13,%xmm14
   0x0000000000001474 <+5236>:	movaps 0x460(%r8),%xmm14
   0x0000000000001488 <+5256>:	movaps 0x4a0(%r8),%xmm15
   0x0000000000001490 <+5264>:	mulps  %xmm10,%xmm14
   0x0000000000001494 <+5268>:	addps  %xmm13,%xmm14
   0x00000000000014c4 <+5316>:	mulps  %xmm10,%xmm15
   0x00000000000014c8 <+5320>:	addps  %xmm14,%xmm15
   0x00000000000014f4 <+5364>:	movaps 0x4e0(%r8),%xmm11
   0x0000000000001500 <+5376>:	mulps  %xmm10,%xmm11
   0x000000000000150c <+5388>:	addps  %xmm7,%xmm11

529
530	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
531	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
532	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000141b <+5147>:	movaps 0x30(%rcx),%xmm12
   0x0000000000001438 <+5176>:	movaps 0x430(%r8),%xmm8
   0x000000000000146c <+5228>:	mulps  %xmm12,%xmm8
   0x0000000000001470 <+5232>:	addps  %xmm14,%xmm8
   0x0000000000001498 <+5272>:	movaps 0x470(%r8),%xmm13
   0x00000000000014a0 <+5280>:	mulps  %xmm12,%xmm13
   0x00000000000014a4 <+5284>:	addps  %xmm14,%xmm13
   0x00000000000014cc <+5324>:	movaps 0x4b0(%r8),%xmm14
   0x00000000000014d4 <+5332>:	mulps  %xmm12,%xmm14
   0x00000000000014d8 <+5336>:	addps  %xmm15,%xmm14
   0x0000000000001510 <+5392>:	movaps 0x4f0(%r8),%xmm7
   0x0000000000001518 <+5400>:	mulps  %xmm12,%xmm7
   0x0000000000001524 <+5412>:	addps  %xmm11,%xmm7

533	      }
534
535	      /* A(2,2)*B(2,1) = C(2,1). */
536	      for (i = 0; i < 4; i++)
537	      {
538	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
539	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
540	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000151c <+5404>:	movaps 0x500(%r8),%xmm12
   0x0000000000001528 <+5416>:	movaps 0x100(%rcx),%xmm9
   0x0000000000001538 <+5432>:	mulps  %xmm9,%xmm12
   0x000000000000153c <+5436>:	addps  %xmm8,%xmm12
   0x0000000000001570 <+5488>:	movaps 0x540(%r8),%xmm15
   0x0000000000001578 <+5496>:	mulps  %xmm9,%xmm15
   0x000000000000157c <+5500>:	addps  %xmm13,%xmm15
   0x00000000000015b0 <+5552>:	movaps 0x580(%r8),%xmm15
   0x00000000000015b8 <+5560>:	mulps  %xmm9,%xmm15
   0x00000000000015bc <+5564>:	addps  %xmm14,%xmm15
   0x00000000000015f0 <+5616>:	movaps 0x5c0(%r8),%xmm15
   0x00000000000015f8 <+5624>:	mulps  %xmm9,%xmm15
   0x0000000000001608 <+5640>:	addps  %xmm7,%xmm15

541
542	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
543	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
544	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001504 <+5380>:	movaps 0x510(%r8),%xmm10
   0x0000000000001530 <+5424>:	movaps 0x110(%rcx),%xmm11
   0x0000000000001548 <+5448>:	mulps  %xmm11,%xmm10
   0x000000000000154c <+5452>:	addps  %xmm12,%xmm10
   0x0000000000001580 <+5504>:	movaps 0x550(%r8),%xmm13
   0x0000000000001588 <+5512>:	mulps  %xmm11,%xmm13
   0x000000000000158c <+5516>:	addps  %xmm15,%xmm13
   0x00000000000015c0 <+5568>:	movaps 0x590(%r8),%xmm14
   0x00000000000015c8 <+5576>:	mulps  %xmm11,%xmm14
   0x00000000000015cc <+5580>:	addps  %xmm15,%xmm14
   0x00000000000015fc <+5628>:	movaps 0x5d0(%r8),%xmm9
   0x0000000000001604 <+5636>:	mulps  %xmm11,%xmm9
   0x0000000000001620 <+5664>:	addps  %xmm15,%xmm9

545
546	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
547	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
548	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000014dc <+5340>:	movaps 0x520(%r8),%xmm15
   0x0000000000001540 <+5440>:	movaps 0x120(%rcx),%xmm8
   0x0000000000001558 <+5464>:	mulps  %xmm8,%xmm15
   0x000000000000155c <+5468>:	addps  %xmm10,%xmm15
   0x0000000000001590 <+5520>:	movaps 0x560(%r8),%xmm15
   0x0000000000001598 <+5528>:	mulps  %xmm8,%xmm15
   0x000000000000159c <+5532>:	addps  %xmm13,%xmm15
   0x00000000000015d0 <+5584>:	movaps 0x5a0(%r8),%xmm15
   0x00000000000015d8 <+5592>:	mulps  %xmm8,%xmm15
   0x00000000000015dc <+5596>:	addps  %xmm14,%xmm15
   0x0000000000001624 <+5668>:	movaps 0x5e0(%r8),%xmm15
   0x000000000000162c <+5676>:	mulps  %xmm8,%xmm15
   0x0000000000001630 <+5680>:	addps  %xmm9,%xmm15

549
550	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
551	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
552	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001550 <+5456>:	movaps 0x130(%rcx),%xmm12
   0x0000000000001560 <+5472>:	movaps 0x530(%r8),%xmm10
   0x0000000000001568 <+5480>:	mulps  %xmm12,%xmm10
   0x000000000000156c <+5484>:	addps  %xmm15,%xmm10
   0x00000000000015a0 <+5536>:	movaps 0x570(%r8),%xmm13
   0x00000000000015a8 <+5544>:	mulps  %xmm12,%xmm13
   0x00000000000015ac <+5548>:	addps  %xmm15,%xmm13
   0x00000000000015e0 <+5600>:	movaps 0x5b0(%r8),%xmm14
   0x00000000000015e8 <+5608>:	mulps  %xmm12,%xmm14
   0x00000000000015ec <+5612>:	addps  %xmm15,%xmm14
   0x000000000000160c <+5644>:	movaps 0x5f0(%r8),%xmm7
   0x0000000000001614 <+5652>:	mulps  %xmm12,%xmm7
   0x0000000000001650 <+5712>:	addps  %xmm15,%xmm7

553	      }
554
555	      /* A(2,3)*B(3,1) = C(2,1). */
556	      for (i = 0; i < 4; i++)
557	      {
558	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
559	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
560	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001618 <+5656>:	movaps 0x600(%r8),%xmm12
   0x0000000000001634 <+5684>:	movaps 0x200(%rcx),%xmm11
   0x000000000000164c <+5708>:	mulps  %xmm11,%xmm12
   0x0000000000001660 <+5728>:	addps  %xmm10,%xmm12
   0x0000000000001690 <+5776>:	movaps 0x640(%r8),%xmm15
   0x0000000000001698 <+5784>:	mulps  %xmm11,%xmm15
   0x000000000000169c <+5788>:	addps  %xmm13,%xmm15
   0x00000000000016d0 <+5840>:	movaps 0x680(%r8),%xmm15
   0x00000000000016d8 <+5848>:	mulps  %xmm11,%xmm15
   0x00000000000016dc <+5852>:	addps  %xmm14,%xmm15
   0x0000000000001710 <+5904>:	movaps 0x6c0(%r8),%xmm15
   0x0000000000001718 <+5912>:	mulps  %xmm11,%xmm15
   0x0000000000001728 <+5928>:	addps  %xmm7,%xmm15

561
562	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
563	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
564	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000163c <+5692>:	movaps 0x210(%rcx),%xmm8
   0x0000000000001664 <+5732>:	movaps 0x610(%r8),%xmm10
   0x000000000000166c <+5740>:	mulps  %xmm8,%xmm10
   0x0000000000001670 <+5744>:	addps  %xmm12,%xmm10
   0x00000000000016a0 <+5792>:	movaps 0x650(%r8),%xmm13
   0x00000000000016a8 <+5800>:	mulps  %xmm8,%xmm13
   0x00000000000016ac <+5804>:	addps  %xmm15,%xmm13
   0x00000000000016e0 <+5856>:	movaps 0x690(%r8),%xmm14
   0x00000000000016e8 <+5864>:	mulps  %xmm8,%xmm14
   0x00000000000016ec <+5868>:	addps  %xmm15,%xmm14
   0x000000000000172c <+5932>:	movaps 0x6d0(%r8),%xmm7
   0x0000000000001734 <+5940>:	mulps  %xmm8,%xmm7
   0x000000000000174c <+5964>:	addps  %xmm15,%xmm7

565
566	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
567	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
568	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001644 <+5700>:	movaps 0x220(%rcx),%xmm9
   0x0000000000001654 <+5716>:	movaps 0x620(%r8),%xmm15
   0x000000000000165c <+5724>:	mulps  %xmm9,%xmm15
   0x000000000000167c <+5756>:	addps  %xmm10,%xmm15
   0x00000000000016b0 <+5808>:	movaps 0x660(%r8),%xmm15
   0x00000000000016b8 <+5816>:	mulps  %xmm9,%xmm15
   0x00000000000016bc <+5820>:	addps  %xmm13,%xmm15
   0x00000000000016f0 <+5872>:	movaps 0x6a0(%r8),%xmm15
   0x00000000000016f8 <+5880>:	mulps  %xmm9,%xmm15
   0x00000000000016fc <+5884>:	addps  %xmm14,%xmm15
   0x000000000000171c <+5916>:	movaps 0x6e0(%r8),%xmm11
   0x0000000000001724 <+5924>:	mulps  %xmm9,%xmm11
   0x0000000000001758 <+5976>:	addps  %xmm7,%xmm11

569
570	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
571	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
572	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001674 <+5748>:	movaps 0x630(%r8),%xmm12
   0x0000000000001680 <+5760>:	movaps 0x230(%rcx),%xmm10
   0x0000000000001688 <+5768>:	mulps  %xmm10,%xmm12
   0x000000000000168c <+5772>:	addps  %xmm15,%xmm12
   0x00000000000016c0 <+5824>:	movaps 0x670(%r8),%xmm13
   0x00000000000016c8 <+5832>:	mulps  %xmm10,%xmm13
   0x00000000000016cc <+5836>:	addps  %xmm15,%xmm13
   0x0000000000001700 <+5888>:	movaps 0x6b0(%r8),%xmm14
   0x0000000000001708 <+5896>:	mulps  %xmm10,%xmm14
   0x000000000000170c <+5900>:	addps  %xmm15,%xmm14
   0x0000000000001738 <+5944>:	movaps 0x6f0(%r8),%xmm8
   0x0000000000001740 <+5952>:	mulps  %xmm10,%xmm8
   0x000000000000176f <+5999>:	addps  %xmm11,%xmm8

573	      }
574
575	      /* A(2,4)*B(4,1) = C(2,1). */
576	      for (i = 0; i < 4; i++)
577	      {
578	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
579	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
580	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001744 <+5956>:	movaps 0x700(%r8),%xmm10
   0x000000000000175c <+5980>:	movaps 0x300(%rcx),%xmm9
   0x000000000000176b <+5995>:	mulps  %xmm9,%xmm10
   0x000000000000177f <+6015>:	addps  %xmm12,%xmm10
   0x00000000000017b3 <+6067>:	movaps 0x740(%r8),%xmm15
   0x00000000000017cb <+6091>:	mulps  %xmm9,%xmm15
   0x00000000000017cf <+6095>:	addps  %xmm13,%xmm15
   0x0000000000001803 <+6147>:	movaps 0x780(%r8),%xmm15
   0x000000000000180b <+6155>:	mulps  %xmm9,%xmm15
   0x000000000000180f <+6159>:	addps  %xmm14,%xmm15
   0x0000000000001843 <+6211>:	movaps 0x7c0(%r8),%xmm15
   0x000000000000184b <+6219>:	mulps  %xmm9,%xmm15
   0x0000000000001857 <+6231>:	addps  %xmm8,%xmm15

581
582	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
583	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
584	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001764 <+5988>:	movaps 0x310(%rcx),%xmm7
   0x0000000000001783 <+6019>:	movaps 0x710(%r8),%xmm12
   0x000000000000178b <+6027>:	mulps  %xmm7,%xmm12
   0x000000000000178f <+6031>:	addps  %xmm10,%xmm12
   0x00000000000017d3 <+6099>:	movaps 0x750(%r8),%xmm13
   0x00000000000017db <+6107>:	mulps  %xmm7,%xmm13
   0x00000000000017df <+6111>:	addps  %xmm15,%xmm13
   0x0000000000001813 <+6163>:	movaps 0x790(%r8),%xmm14
   0x000000000000181b <+6171>:	mulps  %xmm7,%xmm14
   0x000000000000181f <+6175>:	addps  %xmm15,%xmm14
   0x000000000000185b <+6235>:	movaps 0x7d0(%r8),%xmm8
   0x0000000000001863 <+6243>:	mulps  %xmm7,%xmm8
   0x0000000000001873 <+6259>:	addps  %xmm15,%xmm8

585
586	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
587	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
588	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001750 <+5968>:	movaps 0x720(%r8),%xmm15
   0x0000000000001773 <+6003>:	movaps 0x320(%rcx),%xmm11
   0x000000000000177b <+6011>:	mulps  %xmm11,%xmm15
   0x000000000000179b <+6043>:	addps  %xmm12,%xmm15
   0x00000000000017e3 <+6115>:	movaps 0x760(%r8),%xmm15
   0x00000000000017eb <+6123>:	mulps  %xmm11,%xmm15
   0x00000000000017ef <+6127>:	addps  %xmm13,%xmm15
   0x0000000000001823 <+6179>:	movaps 0x7a0(%r8),%xmm15
   0x000000000000182b <+6187>:	mulps  %xmm11,%xmm15
   0x000000000000182f <+6191>:	addps  %xmm14,%xmm15
   0x0000000000001867 <+6247>:	movaps 0x7e0(%r8),%xmm7
   0x000000000000186f <+6255>:	mulps  %xmm11,%xmm7
   0x000000000000187b <+6267>:	addps  %xmm8,%xmm7

589
590	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
591	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
592	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001793 <+6035>:	movaps 0x730(%r8),%xmm10
   0x000000000000179f <+6047>:	movaps 0x330(%rcx),%xmm12
   0x00000000000017a7 <+6055>:	mulps  %xmm12,%xmm10
   0x00000000000017ab <+6059>:	addps  %xmm15,%xmm10
   0x00000000000017f3 <+6131>:	movaps 0x770(%r8),%xmm13
   0x00000000000017fb <+6139>:	mulps  %xmm12,%xmm13
   0x00000000000017ff <+6143>:	addps  %xmm15,%xmm13
   0x0000000000001833 <+6195>:	movaps 0x7b0(%r8),%xmm14
   0x000000000000183b <+6203>:	mulps  %xmm12,%xmm14
   0x000000000000183f <+6207>:	addps  %xmm15,%xmm14
   0x000000000000184f <+6223>:	movaps 0x7f0(%r8),%xmm9
   0x0000000000001877 <+6263>:	mulps  %xmm12,%xmm9
   0x00000000000018a7 <+6311>:	addps  %xmm7,%xmm9

593	      }
594
595	      /* Store C(2,1) block. */
596	      for (i = 0; i < 4; i++)
597	      {
598	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000017af <+6063>:	mulps  %xmm2,%xmm10
   0x000000000000187f <+6271>:	mulps  %xmm2,%xmm13
   0x0000000000001893 <+6291>:	mulps  %xmm2,%xmm14
   0x00000000000018ab <+6315>:	mulps  %xmm2,%xmm9

599	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_21]), C_row[i]);
   0x00000000000017bb <+6075>:	addps  0x100(%r9),%xmm10
   0x0000000000001883 <+6275>:	addps  0x110(%r9),%xmm13
   0x0000000000001897 <+6295>:	addps  0x120(%r9),%xmm14
   0x00000000000018af <+6319>:	addps  0x130(%r9),%xmm9

600	        _mm_store_ps(&C[i*4+C_OFFSET_21], C_row[i]);
   0x00000000000017c3 <+6083>:	movaps %xmm10,0x100(%r9)
   0x000000000000188b <+6283>:	movaps %xmm13,0x110(%r9)
   0x000000000000189f <+6303>:	movaps %xmm14,0x120(%r9)
   0x00000000000018b7 <+6327>:	movaps %xmm9,0x130(%r9)

601	      }
602	    }
603
604	    /* Reset C(2,2) matrix accumulators */
605	    C_row[0] = _mm_setzero_ps();
606	    C_row[1] = _mm_setzero_ps();
607	    C_row[2] = _mm_setzero_ps();
608	    C_row[3] = _mm_setzero_ps();
609
610	    if (norm[4]*norm[17] >= tolerance &&
   0x00000000000018bf <+6335>:	movaps %xmm5,%xmm7
   0x00000000000018c2 <+6338>:	mulss  %xmm0,%xmm7
   0x00000000000018c6 <+6342>:	comiss %xmm1,%xmm7
   0x00000000000018c9 <+6345>:	jb     0x1dc3 <stream_kernel+7619>

611	        norm[5]*norm[21] >= tolerance &&
   0x00000000000018cf <+6351>:	movss  0x2c(%rdx,%rsi,1),%xmm7
   0x00000000000018d5 <+6357>:	mulss  0x6c(%rdx,%rsi,1),%xmm7
   0x00000000000018db <+6363>:	comiss %xmm1,%xmm7
   0x00000000000018de <+6366>:	jb     0x1dc3 <stream_kernel+7619>

612	        norm[6]*norm[25] >= tolerance &&
   0x00000000000018e4 <+6372>:	movss  0x30(%rdx,%rsi,1),%xmm7
   0x00000000000018ea <+6378>:	mulss  0x7c(%rdx,%rsi,1),%xmm7
   0x00000000000018f0 <+6384>:	comiss %xmm1,%xmm7
   0x00000000000018f3 <+6387>:	jb     0x1dc3 <stream_kernel+7619>

613	        norm[7]*norm[29] >= tolerance)
   0x00000000000018f9 <+6393>:	movss  0x34(%rdx,%rsi,1),%xmm7
   0x00000000000018ff <+6399>:	mulss  0x8c(%rdx,%rsi,1),%xmm7
   0x0000000000001908 <+6408>:	comiss %xmm1,%xmm7
   0x000000000000190b <+6411>:	jb     0x1dc3 <stream_kernel+7619>

614	    {
615	      /* A(2,1)*B(1,2) = C(2,2). */
616	      for (i = 0; i < 4; i++)
617	      {
618	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
619	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
620	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001911 <+6417>:	movaps 0x40(%rcx),%xmm7
   0x0000000000001924 <+6436>:	movaps 0x400(%r8),%xmm9
   0x0000000000001944 <+6468>:	movaps 0x440(%r8),%xmm15
   0x000000000000194c <+6476>:	mulps  %xmm7,%xmm9
   0x0000000000001958 <+6488>:	movaps 0x480(%r8),%xmm9
   0x0000000000001980 <+6528>:	mulps  %xmm7,%xmm15
   0x00000000000019b4 <+6580>:	mulps  %xmm7,%xmm9
   0x00000000000019c0 <+6592>:	movaps 0x4c0(%r8),%xmm9
   0x00000000000019e8 <+6632>:	mulps  %xmm7,%xmm9

621
622	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
623	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
624	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001915 <+6421>:	movaps 0x50(%rcx),%xmm11
   0x000000000000192c <+6444>:	movaps 0x410(%r8),%xmm13
   0x0000000000001950 <+6480>:	mulps  %xmm11,%xmm13
   0x0000000000001954 <+6484>:	addps  %xmm9,%xmm13
   0x0000000000001968 <+6504>:	movaps 0x450(%r8),%xmm13
   0x0000000000001984 <+6532>:	mulps  %xmm11,%xmm13
   0x0000000000001988 <+6536>:	addps  %xmm15,%xmm13
   0x00000000000019ac <+6572>:	movaps 0x490(%r8),%xmm14
   0x00000000000019b8 <+6584>:	mulps  %xmm11,%xmm14
   0x00000000000019bc <+6588>:	addps  %xmm9,%xmm14
   0x00000000000019ec <+6636>:	movaps 0x4d0(%r8),%xmm7
   0x00000000000019f4 <+6644>:	mulps  %xmm11,%xmm7
   0x0000000000001a00 <+6656>:	addps  %xmm9,%xmm7

625
626	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
627	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
628	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000191a <+6426>:	movaps 0x60(%rcx),%xmm10
   0x0000000000001934 <+6452>:	movaps 0x420(%r8),%xmm14
   0x0000000000001960 <+6496>:	mulps  %xmm10,%xmm14
   0x0000000000001964 <+6500>:	addps  %xmm13,%xmm14
   0x0000000000001978 <+6520>:	movaps 0x460(%r8),%xmm14
   0x000000000000198c <+6540>:	movaps 0x4a0(%r8),%xmm15
   0x0000000000001994 <+6548>:	mulps  %xmm10,%xmm14
   0x0000000000001998 <+6552>:	addps  %xmm13,%xmm14
   0x00000000000019c8 <+6600>:	mulps  %xmm10,%xmm15
   0x00000000000019cc <+6604>:	addps  %xmm14,%xmm15
   0x00000000000019f8 <+6648>:	movaps 0x4e0(%r8),%xmm11
   0x0000000000001a04 <+6660>:	mulps  %xmm10,%xmm11
   0x0000000000001a10 <+6672>:	addps  %xmm7,%xmm11

629
630	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
631	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
632	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000191f <+6431>:	movaps 0x70(%rcx),%xmm12
   0x000000000000193c <+6460>:	movaps 0x430(%r8),%xmm8
   0x0000000000001970 <+6512>:	mulps  %xmm12,%xmm8
   0x0000000000001974 <+6516>:	addps  %xmm14,%xmm8
   0x000000000000199c <+6556>:	movaps 0x470(%r8),%xmm13
   0x00000000000019a4 <+6564>:	mulps  %xmm12,%xmm13
   0x00000000000019a8 <+6568>:	addps  %xmm14,%xmm13
   0x00000000000019d0 <+6608>:	movaps 0x4b0(%r8),%xmm14
   0x00000000000019d8 <+6616>:	mulps  %xmm12,%xmm14
   0x00000000000019dc <+6620>:	addps  %xmm15,%xmm14
   0x0000000000001a14 <+6676>:	movaps 0x4f0(%r8),%xmm7
   0x0000000000001a1c <+6684>:	mulps  %xmm12,%xmm7
   0x0000000000001a28 <+6696>:	addps  %xmm11,%xmm7

633	      }
634
635	      /* A(2,2)*B(2,2) = C(2,2). */
636	      for (i = 0; i < 4; i++)
637	      {
638	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
639	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
640	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001a20 <+6688>:	movaps 0x500(%r8),%xmm12
   0x0000000000001a2c <+6700>:	movaps 0x140(%rcx),%xmm9
   0x0000000000001a3c <+6716>:	mulps  %xmm9,%xmm12
   0x0000000000001a40 <+6720>:	addps  %xmm8,%xmm12
   0x0000000000001a74 <+6772>:	movaps 0x540(%r8),%xmm15
   0x0000000000001a7c <+6780>:	mulps  %xmm9,%xmm15
   0x0000000000001a80 <+6784>:	addps  %xmm13,%xmm15
   0x0000000000001ab4 <+6836>:	movaps 0x580(%r8),%xmm15
   0x0000000000001abc <+6844>:	mulps  %xmm9,%xmm15
   0x0000000000001ac0 <+6848>:	addps  %xmm14,%xmm15
   0x0000000000001af4 <+6900>:	movaps 0x5c0(%r8),%xmm15
   0x0000000000001afc <+6908>:	mulps  %xmm9,%xmm15
   0x0000000000001b0c <+6924>:	addps  %xmm7,%xmm15

641
642	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
643	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
644	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001a08 <+6664>:	movaps 0x510(%r8),%xmm10
   0x0000000000001a34 <+6708>:	movaps 0x150(%rcx),%xmm11
   0x0000000000001a4c <+6732>:	mulps  %xmm11,%xmm10
   0x0000000000001a50 <+6736>:	addps  %xmm12,%xmm10
   0x0000000000001a84 <+6788>:	movaps 0x550(%r8),%xmm13
   0x0000000000001a8c <+6796>:	mulps  %xmm11,%xmm13
   0x0000000000001a90 <+6800>:	addps  %xmm15,%xmm13
   0x0000000000001ac4 <+6852>:	movaps 0x590(%r8),%xmm14
   0x0000000000001acc <+6860>:	mulps  %xmm11,%xmm14
   0x0000000000001ad0 <+6864>:	addps  %xmm15,%xmm14
   0x0000000000001b00 <+6912>:	movaps 0x5d0(%r8),%xmm9
   0x0000000000001b08 <+6920>:	mulps  %xmm11,%xmm9
   0x0000000000001b24 <+6948>:	addps  %xmm15,%xmm9

645
646	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
647	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
648	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000019e0 <+6624>:	movaps 0x520(%r8),%xmm15
   0x0000000000001a44 <+6724>:	movaps 0x160(%rcx),%xmm8
   0x0000000000001a5c <+6748>:	mulps  %xmm8,%xmm15
   0x0000000000001a60 <+6752>:	addps  %xmm10,%xmm15
   0x0000000000001a94 <+6804>:	movaps 0x560(%r8),%xmm15
   0x0000000000001a9c <+6812>:	mulps  %xmm8,%xmm15
   0x0000000000001aa0 <+6816>:	addps  %xmm13,%xmm15
   0x0000000000001ad4 <+6868>:	movaps 0x5a0(%r8),%xmm15
   0x0000000000001adc <+6876>:	mulps  %xmm8,%xmm15
   0x0000000000001ae0 <+6880>:	addps  %xmm14,%xmm15
   0x0000000000001b28 <+6952>:	movaps 0x5e0(%r8),%xmm15
   0x0000000000001b30 <+6960>:	mulps  %xmm8,%xmm15
   0x0000000000001b34 <+6964>:	addps  %xmm9,%xmm15

649
650	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
651	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
652	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001a54 <+6740>:	movaps 0x170(%rcx),%xmm12
   0x0000000000001a64 <+6756>:	movaps 0x530(%r8),%xmm10
   0x0000000000001a6c <+6764>:	mulps  %xmm12,%xmm10
   0x0000000000001a70 <+6768>:	addps  %xmm15,%xmm10
   0x0000000000001aa4 <+6820>:	movaps 0x570(%r8),%xmm13
   0x0000000000001aac <+6828>:	mulps  %xmm12,%xmm13
   0x0000000000001ab0 <+6832>:	addps  %xmm15,%xmm13
   0x0000000000001ae4 <+6884>:	movaps 0x5b0(%r8),%xmm14
   0x0000000000001aec <+6892>:	mulps  %xmm12,%xmm14
   0x0000000000001af0 <+6896>:	addps  %xmm15,%xmm14
   0x0000000000001b10 <+6928>:	movaps 0x5f0(%r8),%xmm7
   0x0000000000001b18 <+6936>:	mulps  %xmm12,%xmm7
   0x0000000000001b54 <+6996>:	addps  %xmm15,%xmm7

653	      }
654
655	      /* A(2,3)*B(3,2) = C(2,2). */
656	      for (i = 0; i < 4; i++)
657	      {
658	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
659	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
660	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b1c <+6940>:	movaps 0x600(%r8),%xmm12
   0x0000000000001b38 <+6968>:	movaps 0x240(%rcx),%xmm11
   0x0000000000001b50 <+6992>:	mulps  %xmm11,%xmm12
   0x0000000000001b64 <+7012>:	addps  %xmm10,%xmm12
   0x0000000000001b94 <+7060>:	movaps 0x640(%r8),%xmm15
   0x0000000000001b9c <+7068>:	mulps  %xmm11,%xmm15
   0x0000000000001ba0 <+7072>:	addps  %xmm13,%xmm15
   0x0000000000001bd4 <+7124>:	movaps 0x680(%r8),%xmm15
   0x0000000000001bdc <+7132>:	mulps  %xmm11,%xmm15
   0x0000000000001be0 <+7136>:	addps  %xmm14,%xmm15
   0x0000000000001c14 <+7188>:	movaps 0x6c0(%r8),%xmm15
   0x0000000000001c1c <+7196>:	mulps  %xmm11,%xmm15
   0x0000000000001c2c <+7212>:	addps  %xmm7,%xmm15

661
662	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
663	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
664	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b40 <+6976>:	movaps 0x250(%rcx),%xmm8
   0x0000000000001b68 <+7016>:	movaps 0x610(%r8),%xmm10
   0x0000000000001b70 <+7024>:	mulps  %xmm8,%xmm10
   0x0000000000001b74 <+7028>:	addps  %xmm12,%xmm10
   0x0000000000001ba4 <+7076>:	movaps 0x650(%r8),%xmm13
   0x0000000000001bac <+7084>:	mulps  %xmm8,%xmm13
   0x0000000000001bb0 <+7088>:	addps  %xmm15,%xmm13
   0x0000000000001be4 <+7140>:	movaps 0x690(%r8),%xmm14
   0x0000000000001bec <+7148>:	mulps  %xmm8,%xmm14
   0x0000000000001bf0 <+7152>:	addps  %xmm15,%xmm14
   0x0000000000001c30 <+7216>:	movaps 0x6d0(%r8),%xmm7
   0x0000000000001c38 <+7224>:	mulps  %xmm8,%xmm7
   0x0000000000001c50 <+7248>:	addps  %xmm15,%xmm7

665
666	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
667	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
668	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b48 <+6984>:	movaps 0x260(%rcx),%xmm9
   0x0000000000001b58 <+7000>:	movaps 0x620(%r8),%xmm15
   0x0000000000001b60 <+7008>:	mulps  %xmm9,%xmm15
   0x0000000000001b80 <+7040>:	addps  %xmm10,%xmm15
   0x0000000000001bb4 <+7092>:	movaps 0x660(%r8),%xmm15
   0x0000000000001bbc <+7100>:	mulps  %xmm9,%xmm15
   0x0000000000001bc0 <+7104>:	addps  %xmm13,%xmm15
   0x0000000000001bf4 <+7156>:	movaps 0x6a0(%r8),%xmm15
   0x0000000000001bfc <+7164>:	mulps  %xmm9,%xmm15
   0x0000000000001c00 <+7168>:	addps  %xmm14,%xmm15
   0x0000000000001c20 <+7200>:	movaps 0x6e0(%r8),%xmm11
   0x0000000000001c28 <+7208>:	mulps  %xmm9,%xmm11
   0x0000000000001c5c <+7260>:	addps  %xmm7,%xmm11

669
670	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
671	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
672	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b78 <+7032>:	movaps 0x630(%r8),%xmm12
   0x0000000000001b84 <+7044>:	movaps 0x270(%rcx),%xmm10
   0x0000000000001b8c <+7052>:	mulps  %xmm10,%xmm12
   0x0000000000001b90 <+7056>:	addps  %xmm15,%xmm12
   0x0000000000001bc4 <+7108>:	movaps 0x670(%r8),%xmm13
   0x0000000000001bcc <+7116>:	mulps  %xmm10,%xmm13
   0x0000000000001bd0 <+7120>:	addps  %xmm15,%xmm13
   0x0000000000001c04 <+7172>:	movaps 0x6b0(%r8),%xmm14
   0x0000000000001c0c <+7180>:	mulps  %xmm10,%xmm14
   0x0000000000001c10 <+7184>:	addps  %xmm15,%xmm14
   0x0000000000001c3c <+7228>:	movaps 0x6f0(%r8),%xmm8
   0x0000000000001c44 <+7236>:	mulps  %xmm10,%xmm8
   0x0000000000001c73 <+7283>:	addps  %xmm11,%xmm8

673	      }
674
675	      /* A(2,4)*B(4,2) = C(2,2). */
676	      for (i = 0; i < 4; i++)
677	      {
678	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
679	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
680	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c48 <+7240>:	movaps 0x700(%r8),%xmm10
   0x0000000000001c60 <+7264>:	movaps 0x340(%rcx),%xmm9
   0x0000000000001c6f <+7279>:	mulps  %xmm9,%xmm10
   0x0000000000001c83 <+7299>:	addps  %xmm12,%xmm10
   0x0000000000001cb7 <+7351>:	movaps 0x740(%r8),%xmm15
   0x0000000000001ccf <+7375>:	mulps  %xmm9,%xmm15
   0x0000000000001cd3 <+7379>:	addps  %xmm13,%xmm15
   0x0000000000001d07 <+7431>:	movaps 0x780(%r8),%xmm15
   0x0000000000001d0f <+7439>:	mulps  %xmm9,%xmm15
   0x0000000000001d13 <+7443>:	addps  %xmm14,%xmm15
   0x0000000000001d47 <+7495>:	movaps 0x7c0(%r8),%xmm15
   0x0000000000001d4f <+7503>:	mulps  %xmm9,%xmm15
   0x0000000000001d5b <+7515>:	addps  %xmm8,%xmm15

681
682	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
683	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
684	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c68 <+7272>:	movaps 0x350(%rcx),%xmm7
   0x0000000000001c87 <+7303>:	movaps 0x710(%r8),%xmm12
   0x0000000000001c8f <+7311>:	mulps  %xmm7,%xmm12
   0x0000000000001c93 <+7315>:	addps  %xmm10,%xmm12
   0x0000000000001cd7 <+7383>:	movaps 0x750(%r8),%xmm13
   0x0000000000001cdf <+7391>:	mulps  %xmm7,%xmm13
   0x0000000000001ce3 <+7395>:	addps  %xmm15,%xmm13
   0x0000000000001d17 <+7447>:	movaps 0x790(%r8),%xmm14
   0x0000000000001d1f <+7455>:	mulps  %xmm7,%xmm14
   0x0000000000001d23 <+7459>:	addps  %xmm15,%xmm14
   0x0000000000001d5f <+7519>:	movaps 0x7d0(%r8),%xmm8
   0x0000000000001d67 <+7527>:	mulps  %xmm7,%xmm8
   0x0000000000001d77 <+7543>:	addps  %xmm15,%xmm8

685
686	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
687	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
688	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c54 <+7252>:	movaps 0x720(%r8),%xmm15
   0x0000000000001c77 <+7287>:	movaps 0x360(%rcx),%xmm11
   0x0000000000001c7f <+7295>:	mulps  %xmm11,%xmm15
   0x0000000000001c9f <+7327>:	addps  %xmm12,%xmm15
   0x0000000000001ce7 <+7399>:	movaps 0x760(%r8),%xmm15
   0x0000000000001cef <+7407>:	mulps  %xmm11,%xmm15
   0x0000000000001cf3 <+7411>:	addps  %xmm13,%xmm15
   0x0000000000001d27 <+7463>:	movaps 0x7a0(%r8),%xmm15
   0x0000000000001d2f <+7471>:	mulps  %xmm11,%xmm15
   0x0000000000001d33 <+7475>:	addps  %xmm14,%xmm15
   0x0000000000001d6b <+7531>:	movaps 0x7e0(%r8),%xmm7
   0x0000000000001d73 <+7539>:	mulps  %xmm11,%xmm7
   0x0000000000001d7f <+7551>:	addps  %xmm8,%xmm7

689
690	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
691	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
692	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c97 <+7319>:	movaps 0x730(%r8),%xmm10
   0x0000000000001ca3 <+7331>:	movaps 0x370(%rcx),%xmm12
   0x0000000000001cab <+7339>:	mulps  %xmm12,%xmm10
   0x0000000000001caf <+7343>:	addps  %xmm15,%xmm10
   0x0000000000001cf7 <+7415>:	movaps 0x770(%r8),%xmm13
   0x0000000000001cff <+7423>:	mulps  %xmm12,%xmm13
   0x0000000000001d03 <+7427>:	addps  %xmm15,%xmm13
   0x0000000000001d37 <+7479>:	movaps 0x7b0(%r8),%xmm14
   0x0000000000001d3f <+7487>:	mulps  %xmm12,%xmm14
   0x0000000000001d43 <+7491>:	addps  %xmm15,%xmm14
   0x0000000000001d53 <+7507>:	movaps 0x7f0(%r8),%xmm9
   0x0000000000001d7b <+7547>:	mulps  %xmm12,%xmm9
   0x0000000000001dab <+7595>:	addps  %xmm7,%xmm9

693	      }
694
695	      /* Store C(2,2) block. */
696	      for (i = 0; i < 4; i++)
697	      {
698	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000001cb3 <+7347>:	mulps  %xmm2,%xmm10
   0x0000000000001d83 <+7555>:	mulps  %xmm2,%xmm13
   0x0000000000001d97 <+7575>:	mulps  %xmm2,%xmm14
   0x0000000000001daf <+7599>:	mulps  %xmm2,%xmm9

699	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_22]), C_row[i]);
   0x0000000000001cbf <+7359>:	addps  0x140(%r9),%xmm10
   0x0000000000001d87 <+7559>:	addps  0x150(%r9),%xmm13
   0x0000000000001d9b <+7579>:	addps  0x160(%r9),%xmm14
   0x0000000000001db3 <+7603>:	addps  0x170(%r9),%xmm9

700	        _mm_store_ps(&C[i*4+C_OFFSET_22], C_row[i]);
   0x0000000000001cc7 <+7367>:	movaps %xmm10,0x140(%r9)
   0x0000000000001d8f <+7567>:	movaps %xmm13,0x150(%r9)
   0x0000000000001da3 <+7587>:	movaps %xmm14,0x160(%r9)
   0x0000000000001dbb <+7611>:	movaps %xmm9,0x170(%r9)

701	      }
702	    }
703
704	    /* Reset C(2,3) matrix accumulators */
705	    C_row[0] = _mm_setzero_ps();
706	    C_row[1] = _mm_setzero_ps();
707	    C_row[2] = _mm_setzero_ps();
708	    C_row[3] = _mm_setzero_ps();
709
710	    if (norm[4]*norm[18] >= tolerance &&
   0x0000000000001dc3 <+7619>:	movaps %xmm4,%xmm7
   0x0000000000001dc6 <+7622>:	mulss  %xmm0,%xmm7
   0x0000000000001dca <+7626>:	comiss %xmm1,%xmm7
   0x0000000000001dcd <+7629>:	jb     0x22d6 <stream_kernel+8918>

711	        norm[5]*norm[22] >= tolerance &&
   0x0000000000001dd3 <+7635>:	movss  0x2c(%rdx,%rsi,1),%xmm7
   0x0000000000001dd9 <+7641>:	mulss  0x70(%rdx,%rsi,1),%xmm7
   0x0000000000001ddf <+7647>:	comiss %xmm1,%xmm7
   0x0000000000001de2 <+7650>:	jb     0x22d6 <stream_kernel+8918>

712	        norm[6]*norm[26] >= tolerance &&
   0x0000000000001de8 <+7656>:	movss  0x30(%rdx,%rsi,1),%xmm7
   0x0000000000001dee <+7662>:	mulss  0x80(%rdx,%rsi,1),%xmm7
   0x0000000000001df7 <+7671>:	comiss %xmm1,%xmm7
   0x0000000000001dfa <+7674>:	jb     0x22d6 <stream_kernel+8918>

713	        norm[7]*norm[30] >= tolerance)
   0x0000000000001e00 <+7680>:	movss  0x34(%rdx,%rsi,1),%xmm7
   0x0000000000001e06 <+7686>:	mulss  0x90(%rdx,%rsi,1),%xmm7
   0x0000000000001e0f <+7695>:	comiss %xmm1,%xmm7
   0x0000000000001e12 <+7698>:	jb     0x22d6 <stream_kernel+8918>

714	    {
715	      /* A(2,1)*B(1,3) = C(2,3). */
716	      for (i = 0; i < 4; i++)
717	      {
718	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
719	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
720	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001e18 <+7704>:	movaps 0x80(%rcx),%xmm7
   0x0000000000001e37 <+7735>:	movaps 0x400(%r8),%xmm9
   0x0000000000001e57 <+7767>:	movaps 0x440(%r8),%xmm15
   0x0000000000001e5f <+7775>:	mulps  %xmm7,%xmm9
   0x0000000000001e6b <+7787>:	movaps 0x480(%r8),%xmm9
   0x0000000000001e93 <+7827>:	mulps  %xmm7,%xmm15
   0x0000000000001ec7 <+7879>:	mulps  %xmm7,%xmm9
   0x0000000000001ed3 <+7891>:	movaps 0x4c0(%r8),%xmm9
   0x0000000000001efb <+7931>:	mulps  %xmm7,%xmm9

721
722	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
723	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
724	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001e1f <+7711>:	movaps 0x90(%rcx),%xmm11
   0x0000000000001e3f <+7743>:	movaps 0x410(%r8),%xmm13
   0x0000000000001e63 <+7779>:	mulps  %xmm11,%xmm13
   0x0000000000001e67 <+7783>:	addps  %xmm9,%xmm13
   0x0000000000001e7b <+7803>:	movaps 0x450(%r8),%xmm13
   0x0000000000001e97 <+7831>:	mulps  %xmm11,%xmm13
   0x0000000000001e9b <+7835>:	addps  %xmm15,%xmm13
   0x0000000000001ebf <+7871>:	movaps 0x490(%r8),%xmm14
   0x0000000000001ecb <+7883>:	mulps  %xmm11,%xmm14
   0x0000000000001ecf <+7887>:	addps  %xmm9,%xmm14
   0x0000000000001eff <+7935>:	movaps 0x4d0(%r8),%xmm7
   0x0000000000001f07 <+7943>:	mulps  %xmm11,%xmm7
   0x0000000000001f13 <+7955>:	addps  %xmm9,%xmm7

725
726	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
727	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
728	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001e27 <+7719>:	movaps 0xa0(%rcx),%xmm10
   0x0000000000001e47 <+7751>:	movaps 0x420(%r8),%xmm14
   0x0000000000001e73 <+7795>:	mulps  %xmm10,%xmm14
   0x0000000000001e77 <+7799>:	addps  %xmm13,%xmm14
   0x0000000000001e8b <+7819>:	movaps 0x460(%r8),%xmm14
   0x0000000000001e9f <+7839>:	movaps 0x4a0(%r8),%xmm15
   0x0000000000001ea7 <+7847>:	mulps  %xmm10,%xmm14
   0x0000000000001eab <+7851>:	addps  %xmm13,%xmm14
   0x0000000000001edb <+7899>:	mulps  %xmm10,%xmm15
   0x0000000000001edf <+7903>:	addps  %xmm14,%xmm15
   0x0000000000001f0b <+7947>:	movaps 0x4e0(%r8),%xmm11
   0x0000000000001f17 <+7959>:	mulps  %xmm10,%xmm11
   0x0000000000001f23 <+7971>:	addps  %xmm7,%xmm11

729
730	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
731	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
732	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001e2f <+7727>:	movaps 0xb0(%rcx),%xmm12
   0x0000000000001e4f <+7759>:	movaps 0x430(%r8),%xmm8
   0x0000000000001e83 <+7811>:	mulps  %xmm12,%xmm8
   0x0000000000001e87 <+7815>:	addps  %xmm14,%xmm8
   0x0000000000001eaf <+7855>:	movaps 0x470(%r8),%xmm13
   0x0000000000001eb7 <+7863>:	mulps  %xmm12,%xmm13
   0x0000000000001ebb <+7867>:	addps  %xmm14,%xmm13
   0x0000000000001ee3 <+7907>:	movaps 0x4b0(%r8),%xmm14
   0x0000000000001eeb <+7915>:	mulps  %xmm12,%xmm14
   0x0000000000001eef <+7919>:	addps  %xmm15,%xmm14
   0x0000000000001f27 <+7975>:	movaps 0x4f0(%r8),%xmm7
   0x0000000000001f2f <+7983>:	mulps  %xmm12,%xmm7
   0x0000000000001f3b <+7995>:	addps  %xmm11,%xmm7

733	      }
734
735	      /* A(2,2)*B(2,3) = C(2,3). */
736	      for (i = 0; i < 4; i++)
737	      {
738	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
739	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
740	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001f33 <+7987>:	movaps 0x500(%r8),%xmm12
   0x0000000000001f3f <+7999>:	movaps 0x180(%rcx),%xmm9
   0x0000000000001f4f <+8015>:	mulps  %xmm9,%xmm12
   0x0000000000001f53 <+8019>:	addps  %xmm8,%xmm12
   0x0000000000001f87 <+8071>:	movaps 0x540(%r8),%xmm15
   0x0000000000001f8f <+8079>:	mulps  %xmm9,%xmm15
   0x0000000000001f93 <+8083>:	addps  %xmm13,%xmm15
   0x0000000000001fc7 <+8135>:	movaps 0x580(%r8),%xmm15
   0x0000000000001fcf <+8143>:	mulps  %xmm9,%xmm15
   0x0000000000001fd3 <+8147>:	addps  %xmm14,%xmm15
   0x0000000000002007 <+8199>:	movaps 0x5c0(%r8),%xmm15
   0x000000000000200f <+8207>:	mulps  %xmm9,%xmm15
   0x000000000000201f <+8223>:	addps  %xmm7,%xmm15

741
742	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
743	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
744	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001f1b <+7963>:	movaps 0x510(%r8),%xmm10
   0x0000000000001f47 <+8007>:	movaps 0x190(%rcx),%xmm11
   0x0000000000001f5f <+8031>:	mulps  %xmm11,%xmm10
   0x0000000000001f63 <+8035>:	addps  %xmm12,%xmm10
   0x0000000000001f97 <+8087>:	movaps 0x550(%r8),%xmm13
   0x0000000000001f9f <+8095>:	mulps  %xmm11,%xmm13
   0x0000000000001fa3 <+8099>:	addps  %xmm15,%xmm13
   0x0000000000001fd7 <+8151>:	movaps 0x590(%r8),%xmm14
   0x0000000000001fdf <+8159>:	mulps  %xmm11,%xmm14
   0x0000000000001fe3 <+8163>:	addps  %xmm15,%xmm14
   0x0000000000002013 <+8211>:	movaps 0x5d0(%r8),%xmm9
   0x000000000000201b <+8219>:	mulps  %xmm11,%xmm9
   0x0000000000002037 <+8247>:	addps  %xmm15,%xmm9

745
746	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
747	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
748	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001ef3 <+7923>:	movaps 0x520(%r8),%xmm15
   0x0000000000001f57 <+8023>:	movaps 0x1a0(%rcx),%xmm8
   0x0000000000001f6f <+8047>:	mulps  %xmm8,%xmm15
   0x0000000000001f73 <+8051>:	addps  %xmm10,%xmm15
   0x0000000000001fa7 <+8103>:	movaps 0x560(%r8),%xmm15
   0x0000000000001faf <+8111>:	mulps  %xmm8,%xmm15
   0x0000000000001fb3 <+8115>:	addps  %xmm13,%xmm15
   0x0000000000001fe7 <+8167>:	movaps 0x5a0(%r8),%xmm15
   0x0000000000001fef <+8175>:	mulps  %xmm8,%xmm15
   0x0000000000001ff3 <+8179>:	addps  %xmm14,%xmm15
   0x000000000000203b <+8251>:	movaps 0x5e0(%r8),%xmm15
   0x0000000000002043 <+8259>:	mulps  %xmm8,%xmm15
   0x0000000000002047 <+8263>:	addps  %xmm9,%xmm15

749
750	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
751	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
752	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001f67 <+8039>:	movaps 0x1b0(%rcx),%xmm12
   0x0000000000001f77 <+8055>:	movaps 0x530(%r8),%xmm10
   0x0000000000001f7f <+8063>:	mulps  %xmm12,%xmm10
   0x0000000000001f83 <+8067>:	addps  %xmm15,%xmm10
   0x0000000000001fb7 <+8119>:	movaps 0x570(%r8),%xmm13
   0x0000000000001fbf <+8127>:	mulps  %xmm12,%xmm13
   0x0000000000001fc3 <+8131>:	addps  %xmm15,%xmm13
   0x0000000000001ff7 <+8183>:	movaps 0x5b0(%r8),%xmm14
   0x0000000000001fff <+8191>:	mulps  %xmm12,%xmm14
   0x0000000000002003 <+8195>:	addps  %xmm15,%xmm14
   0x0000000000002023 <+8227>:	movaps 0x5f0(%r8),%xmm7
   0x000000000000202b <+8235>:	mulps  %xmm12,%xmm7
   0x0000000000002067 <+8295>:	addps  %xmm15,%xmm7

753	      }
754
755	      /* A(2,3)*B(3,3) = C(2,3). */
756	      for (i = 0; i < 4; i++)
757	      {
758	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
759	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
760	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000202f <+8239>:	movaps 0x600(%r8),%xmm12
   0x000000000000204b <+8267>:	movaps 0x280(%rcx),%xmm11
   0x0000000000002063 <+8291>:	mulps  %xmm11,%xmm12
   0x0000000000002077 <+8311>:	addps  %xmm10,%xmm12
   0x00000000000020a7 <+8359>:	movaps 0x640(%r8),%xmm15
   0x00000000000020af <+8367>:	mulps  %xmm11,%xmm15
   0x00000000000020b3 <+8371>:	addps  %xmm13,%xmm15
   0x00000000000020e7 <+8423>:	movaps 0x680(%r8),%xmm15
   0x00000000000020ef <+8431>:	mulps  %xmm11,%xmm15
   0x00000000000020f3 <+8435>:	addps  %xmm14,%xmm15
   0x0000000000002127 <+8487>:	movaps 0x6c0(%r8),%xmm15
   0x000000000000212f <+8495>:	mulps  %xmm11,%xmm15
   0x000000000000213f <+8511>:	addps  %xmm7,%xmm15

761
762	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
763	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
764	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002053 <+8275>:	movaps 0x290(%rcx),%xmm8
   0x000000000000207b <+8315>:	movaps 0x610(%r8),%xmm10
   0x0000000000002083 <+8323>:	mulps  %xmm8,%xmm10
   0x0000000000002087 <+8327>:	addps  %xmm12,%xmm10
   0x00000000000020b7 <+8375>:	movaps 0x650(%r8),%xmm13
   0x00000000000020bf <+8383>:	mulps  %xmm8,%xmm13
   0x00000000000020c3 <+8387>:	addps  %xmm15,%xmm13
   0x00000000000020f7 <+8439>:	movaps 0x690(%r8),%xmm14
   0x00000000000020ff <+8447>:	mulps  %xmm8,%xmm14
   0x0000000000002103 <+8451>:	addps  %xmm15,%xmm14
   0x0000000000002143 <+8515>:	movaps 0x6d0(%r8),%xmm7
   0x000000000000214b <+8523>:	mulps  %xmm8,%xmm7
   0x0000000000002163 <+8547>:	addps  %xmm15,%xmm7

765
766	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
767	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
768	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000205b <+8283>:	movaps 0x2a0(%rcx),%xmm9
   0x000000000000206b <+8299>:	movaps 0x620(%r8),%xmm15
   0x0000000000002073 <+8307>:	mulps  %xmm9,%xmm15
   0x0000000000002093 <+8339>:	addps  %xmm10,%xmm15
   0x00000000000020c7 <+8391>:	movaps 0x660(%r8),%xmm15
   0x00000000000020cf <+8399>:	mulps  %xmm9,%xmm15
   0x00000000000020d3 <+8403>:	addps  %xmm13,%xmm15
   0x0000000000002107 <+8455>:	movaps 0x6a0(%r8),%xmm15
   0x000000000000210f <+8463>:	mulps  %xmm9,%xmm15
   0x0000000000002113 <+8467>:	addps  %xmm14,%xmm15
   0x0000000000002133 <+8499>:	movaps 0x6e0(%r8),%xmm11
   0x000000000000213b <+8507>:	mulps  %xmm9,%xmm11
   0x000000000000216f <+8559>:	addps  %xmm7,%xmm11

769
770	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
771	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
772	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000208b <+8331>:	movaps 0x630(%r8),%xmm12
   0x0000000000002097 <+8343>:	movaps 0x2b0(%rcx),%xmm10
   0x000000000000209f <+8351>:	mulps  %xmm10,%xmm12
   0x00000000000020a3 <+8355>:	addps  %xmm15,%xmm12
   0x00000000000020d7 <+8407>:	movaps 0x670(%r8),%xmm13
   0x00000000000020df <+8415>:	mulps  %xmm10,%xmm13
   0x00000000000020e3 <+8419>:	addps  %xmm15,%xmm13
   0x0000000000002117 <+8471>:	movaps 0x6b0(%r8),%xmm14
   0x000000000000211f <+8479>:	mulps  %xmm10,%xmm14
   0x0000000000002123 <+8483>:	addps  %xmm15,%xmm14
   0x000000000000214f <+8527>:	movaps 0x6f0(%r8),%xmm8
   0x0000000000002157 <+8535>:	mulps  %xmm10,%xmm8
   0x0000000000002186 <+8582>:	addps  %xmm11,%xmm8

773	      }
774
775	      /* A(2,4)*B(4,3) = C(2,3). */
776	      for (i = 0; i < 4; i++)
777	      {
778	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
779	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
780	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000215b <+8539>:	movaps 0x700(%r8),%xmm10
   0x0000000000002173 <+8563>:	movaps 0x380(%rcx),%xmm9
   0x0000000000002182 <+8578>:	mulps  %xmm9,%xmm10
   0x0000000000002196 <+8598>:	addps  %xmm12,%xmm10
   0x00000000000021ca <+8650>:	movaps 0x740(%r8),%xmm15
   0x00000000000021e2 <+8674>:	mulps  %xmm9,%xmm15
   0x00000000000021e6 <+8678>:	addps  %xmm13,%xmm15
   0x000000000000221a <+8730>:	movaps 0x780(%r8),%xmm15
   0x0000000000002222 <+8738>:	mulps  %xmm9,%xmm15
   0x0000000000002226 <+8742>:	addps  %xmm14,%xmm15
   0x000000000000225a <+8794>:	movaps 0x7c0(%r8),%xmm15
   0x0000000000002262 <+8802>:	mulps  %xmm9,%xmm15
   0x000000000000226e <+8814>:	addps  %xmm8,%xmm15

781
782	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
783	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
784	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000217b <+8571>:	movaps 0x390(%rcx),%xmm7
   0x000000000000219a <+8602>:	movaps 0x710(%r8),%xmm12
   0x00000000000021a2 <+8610>:	mulps  %xmm7,%xmm12
   0x00000000000021a6 <+8614>:	addps  %xmm10,%xmm12
   0x00000000000021ea <+8682>:	movaps 0x750(%r8),%xmm13
   0x00000000000021f2 <+8690>:	mulps  %xmm7,%xmm13
   0x00000000000021f6 <+8694>:	addps  %xmm15,%xmm13
   0x000000000000222a <+8746>:	movaps 0x790(%r8),%xmm14
   0x0000000000002232 <+8754>:	mulps  %xmm7,%xmm14
   0x0000000000002236 <+8758>:	addps  %xmm15,%xmm14
   0x0000000000002272 <+8818>:	movaps 0x7d0(%r8),%xmm8
   0x000000000000227a <+8826>:	mulps  %xmm7,%xmm8
   0x000000000000228a <+8842>:	addps  %xmm15,%xmm8

785
786	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
787	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
788	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002167 <+8551>:	movaps 0x720(%r8),%xmm15
   0x000000000000218a <+8586>:	movaps 0x3a0(%rcx),%xmm11
   0x0000000000002192 <+8594>:	mulps  %xmm11,%xmm15
   0x00000000000021b2 <+8626>:	addps  %xmm12,%xmm15
   0x00000000000021fa <+8698>:	movaps 0x760(%r8),%xmm15
   0x0000000000002202 <+8706>:	mulps  %xmm11,%xmm15
   0x0000000000002206 <+8710>:	addps  %xmm13,%xmm15
   0x000000000000223a <+8762>:	movaps 0x7a0(%r8),%xmm15
   0x0000000000002242 <+8770>:	mulps  %xmm11,%xmm15
   0x0000000000002246 <+8774>:	addps  %xmm14,%xmm15
   0x000000000000227e <+8830>:	movaps 0x7e0(%r8),%xmm7
   0x0000000000002286 <+8838>:	mulps  %xmm11,%xmm7
   0x0000000000002292 <+8850>:	addps  %xmm8,%xmm7

789
790	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
791	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
792	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000021aa <+8618>:	movaps 0x730(%r8),%xmm10
   0x00000000000021b6 <+8630>:	movaps 0x3b0(%rcx),%xmm12
   0x00000000000021be <+8638>:	mulps  %xmm12,%xmm10
   0x00000000000021c2 <+8642>:	addps  %xmm15,%xmm10
   0x000000000000220a <+8714>:	movaps 0x770(%r8),%xmm13
   0x0000000000002212 <+8722>:	mulps  %xmm12,%xmm13
   0x0000000000002216 <+8726>:	addps  %xmm15,%xmm13
   0x000000000000224a <+8778>:	movaps 0x7b0(%r8),%xmm14
   0x0000000000002252 <+8786>:	mulps  %xmm12,%xmm14
   0x0000000000002256 <+8790>:	addps  %xmm15,%xmm14
   0x0000000000002266 <+8806>:	movaps 0x7f0(%r8),%xmm9
   0x000000000000228e <+8846>:	mulps  %xmm12,%xmm9
   0x00000000000022be <+8894>:	addps  %xmm7,%xmm9

793	      }
794
795	      /* Store C(2,3) block. */
796	      for (i = 0; i < 4; i++)
797	      {
798	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000021c6 <+8646>:	mulps  %xmm2,%xmm10
   0x0000000000002296 <+8854>:	mulps  %xmm2,%xmm13
   0x00000000000022aa <+8874>:	mulps  %xmm2,%xmm14
   0x00000000000022c2 <+8898>:	mulps  %xmm2,%xmm9

799	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_23]), C_row[i]);
   0x00000000000021d2 <+8658>:	addps  0x180(%r9),%xmm10
   0x000000000000229a <+8858>:	addps  0x190(%r9),%xmm13
   0x00000000000022ae <+8878>:	addps  0x1a0(%r9),%xmm14
   0x00000000000022c6 <+8902>:	addps  0x1b0(%r9),%xmm9

800	        _mm_store_ps(&C[i*4+C_OFFSET_23], C_row[i]);
   0x00000000000021da <+8666>:	movaps %xmm10,0x180(%r9)
   0x00000000000022a2 <+8866>:	movaps %xmm13,0x190(%r9)
   0x00000000000022b6 <+8886>:	movaps %xmm14,0x1a0(%r9)
   0x00000000000022ce <+8910>:	movaps %xmm9,0x1b0(%r9)

801	      }
802	    }
803
804	    /* Reset C(2,4) matrix accumulators */
805	    C_row[0] = _mm_setzero_ps();
806	    C_row[1] = _mm_setzero_ps();
807	    C_row[2] = _mm_setzero_ps();
808	    C_row[3] = _mm_setzero_ps();
809
810	    if (norm[4]*norm[19] >= tolerance &&
   0x00000000000022d6 <+8918>:	mulss  %xmm3,%xmm0
   0x00000000000022da <+8922>:	comiss %xmm1,%xmm0
   0x00000000000022dd <+8925>:	jb     0x27df <stream_kernel+10207>

811	        norm[5]*norm[23] >= tolerance &&
   0x00000000000022e3 <+8931>:	movss  0x2c(%rdx,%rsi,1),%xmm0
   0x00000000000022e9 <+8937>:	mulss  0x74(%rdx,%rsi,1),%xmm0
   0x00000000000022ef <+8943>:	comiss %xmm1,%xmm0
   0x00000000000022f2 <+8946>:	jb     0x27df <stream_kernel+10207>

812	        norm[6]*norm[27] >= tolerance &&
   0x00000000000022f8 <+8952>:	movss  0x30(%rdx,%rsi,1),%xmm0
   0x00000000000022fe <+8958>:	mulss  0x84(%rdx,%rsi,1),%xmm0
   0x0000000000002307 <+8967>:	comiss %xmm1,%xmm0
   0x000000000000230a <+8970>:	jb     0x27df <stream_kernel+10207>

813	        norm[7]*norm[31] >= tolerance)
   0x0000000000002310 <+8976>:	movss  0x34(%rdx,%rsi,1),%xmm0
   0x0000000000002316 <+8982>:	mulss  0x94(%rdx,%rsi,1),%xmm0
   0x000000000000231f <+8991>:	comiss %xmm1,%xmm0
   0x0000000000002322 <+8994>:	jb     0x27df <stream_kernel+10207>

814	    {
815	      /* A(2,1)*B(1,4) = C(2,4). */
816	      for (i = 0; i < 4; i++)
817	      {
818	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
819	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
820	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002328 <+9000>:	movaps 0xc0(%rcx),%xmm7
   0x0000000000002346 <+9030>:	movaps 0x400(%r8),%xmm10
   0x0000000000002366 <+9062>:	movaps 0x440(%r8),%xmm15
   0x0000000000002376 <+9078>:	mulps  %xmm7,%xmm10
   0x00000000000023a2 <+9122>:	movaps 0x480(%r8),%xmm11
   0x00000000000023aa <+9130>:	mulps  %xmm7,%xmm15
   0x00000000000023c6 <+9158>:	movaps 0x4c0(%r8),%xmm14
   0x00000000000023de <+9182>:	mulps  %xmm7,%xmm11
   0x0000000000002412 <+9234>:	mulps  %xmm7,%xmm14

821
822	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
823	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
824	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000232f <+9007>:	movaps 0xd0(%rcx),%xmm0
   0x000000000000234e <+9038>:	movaps 0x410(%r8),%xmm9
   0x000000000000236e <+9070>:	movaps 0x450(%r8),%xmm14
   0x000000000000237a <+9082>:	mulps  %xmm0,%xmm9
   0x000000000000237e <+9086>:	addps  %xmm10,%xmm9
   0x00000000000023ae <+9134>:	mulps  %xmm0,%xmm14
   0x00000000000023b2 <+9138>:	addps  %xmm15,%xmm14
   0x00000000000023b6 <+9142>:	movaps 0x490(%r8),%xmm15
   0x00000000000023e2 <+9186>:	mulps  %xmm0,%xmm15
   0x00000000000023e6 <+9190>:	addps  %xmm11,%xmm15
   0x0000000000002416 <+9238>:	movaps 0x4d0(%r8),%xmm7
   0x000000000000241e <+9246>:	mulps  %xmm0,%xmm7
   0x0000000000002429 <+9257>:	addps  %xmm14,%xmm7

825
826	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
827	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
828	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002336 <+9014>:	movaps 0xe0(%rcx),%xmm8
   0x0000000000002356 <+9046>:	movaps 0x420(%r8),%xmm11
   0x0000000000002382 <+9090>:	movaps 0x460(%r8),%xmm10
   0x000000000000238a <+9098>:	mulps  %xmm8,%xmm11
   0x000000000000238e <+9102>:	addps  %xmm9,%xmm11
   0x00000000000023be <+9150>:	mulps  %xmm8,%xmm10
   0x00000000000023c2 <+9154>:	addps  %xmm14,%xmm10
   0x00000000000023ea <+9194>:	movaps 0x4a0(%r8),%xmm11
   0x00000000000023f2 <+9202>:	mulps  %xmm8,%xmm11
   0x00000000000023f6 <+9206>:	addps  %xmm15,%xmm11
   0x0000000000002421 <+9249>:	movaps 0x4e0(%r8),%xmm0
   0x000000000000242d <+9261>:	mulps  %xmm8,%xmm0
   0x0000000000002445 <+9285>:	addps  %xmm7,%xmm0

829
830	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
831	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
832	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000233e <+9022>:	movaps 0xf0(%rcx),%xmm12
   0x000000000000235e <+9054>:	movaps 0x430(%r8),%xmm13
   0x0000000000002392 <+9106>:	movaps 0x470(%r8),%xmm9
   0x000000000000239a <+9114>:	mulps  %xmm12,%xmm13
   0x000000000000239e <+9118>:	addps  %xmm11,%xmm13
   0x00000000000023ce <+9166>:	mulps  %xmm12,%xmm9
   0x00000000000023d2 <+9170>:	addps  %xmm10,%xmm9
   0x00000000000023d6 <+9174>:	movaps 0x4b0(%r8),%xmm10
   0x0000000000002402 <+9218>:	mulps  %xmm12,%xmm10
   0x0000000000002406 <+9222>:	addps  %xmm11,%xmm10
   0x000000000000240a <+9226>:	movaps 0x4f0(%r8),%xmm11
   0x0000000000002441 <+9281>:	mulps  %xmm12,%xmm11
   0x000000000000245c <+9308>:	addps  %xmm0,%xmm11

833	      }
834
835	      /* A(2,2)*B(2,4) = C(2,4). */
836	      for (i = 0; i < 4; i++)
837	      {
838	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
839	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
840	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002431 <+9265>:	movaps 0x1c0(%rcx),%xmm8
   0x0000000000002448 <+9288>:	movaps 0x500(%r8),%xmm12
   0x0000000000002458 <+9304>:	mulps  %xmm8,%xmm12
   0x0000000000002460 <+9312>:	movaps 0x540(%r8),%xmm0
   0x000000000000246c <+9324>:	addps  %xmm13,%xmm12
   0x0000000000002480 <+9344>:	mulps  %xmm8,%xmm0
   0x000000000000249c <+9372>:	addps  %xmm9,%xmm0
   0x00000000000024d0 <+9424>:	movaps 0x580(%r8),%xmm9
   0x00000000000024d8 <+9432>:	mulps  %xmm8,%xmm9
   0x00000000000024eb <+9451>:	addps  %xmm10,%xmm9
   0x000000000000250f <+9487>:	movaps 0x5c0(%r8),%xmm10
   0x0000000000002517 <+9495>:	mulps  %xmm8,%xmm10
   0x000000000000253b <+9531>:	addps  %xmm11,%xmm10

841
842	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
843	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
844	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002439 <+9273>:	movaps 0x1d0(%rcx),%xmm14
   0x0000000000002470 <+9328>:	movaps 0x510(%r8),%xmm13
   0x0000000000002478 <+9336>:	mulps  %xmm14,%xmm13
   0x000000000000247c <+9340>:	addps  %xmm12,%xmm13
   0x00000000000024a0 <+9376>:	movaps 0x550(%r8),%xmm9
   0x00000000000024a8 <+9384>:	mulps  %xmm14,%xmm9
   0x00000000000024bc <+9404>:	addps  %xmm0,%xmm9
   0x00000000000024ef <+9455>:	movaps 0x590(%r8),%xmm10
   0x00000000000024f7 <+9463>:	mulps  %xmm14,%xmm10
   0x00000000000024fb <+9467>:	addps  %xmm9,%xmm10
   0x000000000000251b <+9499>:	movaps 0x5d0(%r8),%xmm8
   0x0000000000002523 <+9507>:	mulps  %xmm14,%xmm8
   0x0000000000002557 <+9559>:	addps  %xmm10,%xmm8

845
846	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
847	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
848	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000023fa <+9210>:	movaps 0x1e0(%rcx),%xmm15
   0x0000000000002450 <+9296>:	movaps 0x520(%r8),%xmm7
   0x0000000000002468 <+9320>:	mulps  %xmm15,%xmm7
   0x000000000000248c <+9356>:	addps  %xmm13,%xmm7
   0x00000000000024c0 <+9408>:	movaps 0x560(%r8),%xmm0
   0x00000000000024c8 <+9416>:	mulps  %xmm15,%xmm0
   0x00000000000024cc <+9420>:	addps  %xmm9,%xmm0
   0x00000000000024ff <+9471>:	movaps 0x5a0(%r8),%xmm9
   0x0000000000002507 <+9479>:	mulps  %xmm15,%xmm9
   0x000000000000250b <+9483>:	addps  %xmm10,%xmm9
   0x000000000000255b <+9563>:	movaps 0x5e0(%r8),%xmm10
   0x0000000000002563 <+9571>:	mulps  %xmm15,%xmm10
   0x0000000000002577 <+9591>:	addps  %xmm8,%xmm10

849
850	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
851	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
852	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002484 <+9348>:	movaps 0x1f0(%rcx),%xmm12
   0x0000000000002490 <+9360>:	movaps 0x530(%r8),%xmm13
   0x0000000000002498 <+9368>:	mulps  %xmm12,%xmm13
   0x00000000000024ac <+9388>:	addps  %xmm7,%xmm13
   0x00000000000024b0 <+9392>:	movaps 0x570(%r8),%xmm7
   0x00000000000024b8 <+9400>:	mulps  %xmm12,%xmm7
   0x00000000000024dc <+9436>:	addps  %xmm0,%xmm7
   0x00000000000024df <+9439>:	movaps 0x5b0(%r8),%xmm0
   0x00000000000024e7 <+9447>:	mulps  %xmm12,%xmm0
   0x0000000000002527 <+9511>:	addps  %xmm9,%xmm0
   0x000000000000253f <+9535>:	movaps 0x5f0(%r8),%xmm11
   0x0000000000002547 <+9543>:	mulps  %xmm12,%xmm11
   0x0000000000002583 <+9603>:	addps  %xmm10,%xmm11

853	      }
854
855	      /* A(2,3)*B(3,4) = C(2,4). */
856	      for (i = 0; i < 4; i++)
857	      {
858	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
859	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
860	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000252b <+9515>:	movaps 0x2c0(%rcx),%xmm14
   0x000000000000254b <+9547>:	movaps 0x600(%r8),%xmm12
   0x0000000000002553 <+9555>:	mulps  %xmm14,%xmm12
   0x0000000000002567 <+9575>:	addps  %xmm13,%xmm12
   0x000000000000259f <+9631>:	movaps 0x640(%r8),%xmm12
   0x00000000000025a7 <+9639>:	mulps  %xmm14,%xmm12
   0x00000000000025bb <+9659>:	addps  %xmm7,%xmm12
   0x00000000000025df <+9695>:	movaps 0x680(%r8),%xmm12
   0x00000000000025e7 <+9703>:	mulps  %xmm14,%xmm12
   0x00000000000025fb <+9723>:	addps  %xmm0,%xmm12
   0x000000000000261f <+9759>:	movaps 0x6c0(%r8),%xmm12
   0x0000000000002627 <+9767>:	mulps  %xmm14,%xmm12
   0x0000000000002647 <+9799>:	addps  %xmm11,%xmm12

861
862	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
863	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
864	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002533 <+9523>:	movaps 0x2d0(%rcx),%xmm9
   0x000000000000256b <+9579>:	movaps 0x610(%r8),%xmm13
   0x0000000000002573 <+9587>:	mulps  %xmm9,%xmm13
   0x000000000000259b <+9627>:	addps  %xmm12,%xmm13
   0x00000000000025bf <+9663>:	movaps 0x650(%r8),%xmm7
   0x00000000000025c7 <+9671>:	mulps  %xmm9,%xmm7
   0x00000000000025db <+9691>:	addps  %xmm12,%xmm7
   0x00000000000025ff <+9727>:	movaps 0x690(%r8),%xmm0
   0x0000000000002607 <+9735>:	mulps  %xmm9,%xmm0
   0x000000000000261b <+9755>:	addps  %xmm12,%xmm0
   0x000000000000264b <+9803>:	movaps 0x6d0(%r8),%xmm11
   0x0000000000002653 <+9811>:	mulps  %xmm9,%xmm11
   0x0000000000002677 <+9847>:	addps  %xmm12,%xmm11

865
866	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
867	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
868	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000257b <+9595>:	movaps 0x620(%r8),%xmm8
   0x0000000000002587 <+9607>:	movaps 0x2e0(%rcx),%xmm10
   0x0000000000002597 <+9623>:	mulps  %xmm10,%xmm8
   0x00000000000025ab <+9643>:	addps  %xmm13,%xmm8
   0x00000000000025cf <+9679>:	movaps 0x660(%r8),%xmm8
   0x00000000000025d7 <+9687>:	mulps  %xmm10,%xmm8
   0x00000000000025eb <+9707>:	addps  %xmm7,%xmm8
   0x000000000000260f <+9743>:	movaps 0x6a0(%r8),%xmm8
   0x0000000000002617 <+9751>:	mulps  %xmm10,%xmm8
   0x0000000000002637 <+9783>:	addps  %xmm0,%xmm8
   0x0000000000002657 <+9815>:	movaps 0x6e0(%r8),%xmm9
   0x000000000000265f <+9823>:	mulps  %xmm10,%xmm9
   0x0000000000002683 <+9859>:	addps  %xmm11,%xmm9

869
870	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
871	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
872	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000258f <+9615>:	movaps 0x2f0(%rcx),%xmm15
   0x00000000000025af <+9647>:	movaps 0x630(%r8),%xmm13
   0x00000000000025b7 <+9655>:	mulps  %xmm15,%xmm13
   0x00000000000025cb <+9675>:	addps  %xmm8,%xmm13
   0x00000000000025ef <+9711>:	movaps 0x670(%r8),%xmm7
   0x00000000000025f7 <+9719>:	mulps  %xmm15,%xmm7
   0x000000000000260b <+9739>:	addps  %xmm8,%xmm7
   0x000000000000262b <+9771>:	movaps 0x6f0(%r8),%xmm14
   0x0000000000002633 <+9779>:	mulps  %xmm15,%xmm14
   0x000000000000263b <+9787>:	movaps 0x6b0(%r8),%xmm0
   0x0000000000002643 <+9795>:	mulps  %xmm15,%xmm0
   0x000000000000266b <+9835>:	addps  %xmm8,%xmm0
   0x000000000000269b <+9883>:	addps  %xmm9,%xmm14

873	      }
874
875	      /* A(2,4)*B(4,4) = C(2,4). */
876	      for (i = 0; i < 4; i++)
877	      {
878	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
879	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
880	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002663 <+9827>:	movaps 0x700(%r8),%xmm10
   0x0000000000002687 <+9863>:	movaps 0x3c0(%rcx),%xmm15
   0x0000000000002697 <+9879>:	mulps  %xmm15,%xmm10
   0x00000000000026ab <+9899>:	addps  %xmm13,%xmm10
   0x00000000000026cf <+9935>:	movaps 0x740(%r8),%xmm8
   0x00000000000026d7 <+9943>:	mulps  %xmm15,%xmm8
   0x00000000000026eb <+9963>:	addps  %xmm7,%xmm8
   0x000000000000270f <+9999>:	movaps 0x780(%r8),%xmm7
   0x0000000000002717 <+10007>:	mulps  %xmm15,%xmm7
   0x000000000000272b <+10027>:	addps  %xmm0,%xmm7
   0x000000000000274c <+10060>:	movaps 0x7c0(%r8),%xmm0
   0x0000000000002754 <+10068>:	mulps  %xmm15,%xmm0
   0x00000000000027a4 <+10148>:	addps  %xmm14,%xmm0

881
882	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
883	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
884	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000266f <+9839>:	movaps 0x710(%r8),%xmm8
   0x000000000000268f <+9871>:	movaps 0x3d0(%rcx),%xmm11
   0x000000000000269f <+9887>:	mulps  %xmm11,%xmm8
   0x00000000000026bb <+9915>:	addps  %xmm10,%xmm8
   0x00000000000026ef <+9967>:	movaps 0x750(%r8),%xmm7
   0x00000000000026f7 <+9975>:	mulps  %xmm11,%xmm7
   0x00000000000026fb <+9979>:	addps  %xmm8,%xmm7
   0x000000000000272e <+10030>:	movaps 0x790(%r8),%xmm0
   0x0000000000002736 <+10038>:	mulps  %xmm11,%xmm0
   0x000000000000273a <+10042>:	addps  %xmm7,%xmm0
   0x00000000000027a8 <+10152>:	movaps 0x7d0(%r8),%xmm14
   0x00000000000027b0 <+10160>:	mulps  %xmm11,%xmm14
   0x00000000000027b4 <+10164>:	addps  %xmm0,%xmm14

885
886	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
887	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
888	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000267b <+9851>:	movaps 0x720(%r8),%xmm12
   0x00000000000026af <+9903>:	movaps 0x3e0(%rcx),%xmm13
   0x00000000000026b7 <+9911>:	mulps  %xmm13,%xmm12
   0x00000000000026cb <+9931>:	addps  %xmm8,%xmm12
   0x00000000000026df <+9951>:	movaps 0x760(%r8),%xmm12
   0x00000000000026e7 <+9959>:	mulps  %xmm13,%xmm12
   0x000000000000270b <+9995>:	addps  %xmm7,%xmm12
   0x000000000000273d <+10045>:	movaps 0x7a0(%r8),%xmm7
   0x0000000000002745 <+10053>:	mulps  %xmm13,%xmm7
   0x0000000000002749 <+10057>:	addps  %xmm0,%xmm7
   0x0000000000002758 <+10072>:	movaps 0x7e0(%r8),%xmm15
   0x0000000000002760 <+10080>:	mulps  %xmm13,%xmm15
   0x00000000000027c4 <+10180>:	addps  %xmm14,%xmm15

889
890	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
891	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
892	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000026a3 <+9891>:	movaps 0x730(%r8),%xmm9
   0x00000000000026bf <+9919>:	movaps 0x3f0(%rcx),%xmm10
   0x00000000000026c7 <+9927>:	mulps  %xmm10,%xmm9
   0x00000000000026db <+9947>:	addps  %xmm12,%xmm9
   0x00000000000026ff <+9983>:	movaps 0x770(%r8),%xmm8
   0x0000000000002707 <+9991>:	mulps  %xmm10,%xmm8
   0x000000000000271b <+10011>:	addps  %xmm12,%xmm8
   0x000000000000271f <+10015>:	movaps 0x7b0(%r8),%xmm12
   0x0000000000002727 <+10023>:	mulps  %xmm10,%xmm12
   0x0000000000002764 <+10084>:	addps  %xmm7,%xmm12
   0x00000000000027b8 <+10168>:	movaps 0x7f0(%r8),%xmm0
   0x00000000000027c0 <+10176>:	mulps  %xmm10,%xmm0
   0x00000000000027c8 <+10184>:	addps  %xmm15,%xmm0

893	      }
894
895	      /* Store C(2,4) block. */
896	      for (i = 0; i < 4; i++)
897	      {
898	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000002768 <+10088>:	mulps  %xmm2,%xmm9
   0x000000000000277c <+10108>:	mulps  %xmm2,%xmm8
   0x0000000000002790 <+10128>:	mulps  %xmm2,%xmm12
   0x00000000000027cc <+10188>:	mulps  %xmm2,%xmm0

899	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_24]), C_row[i]);
   0x000000000000276c <+10092>:	addps  0x1c0(%r9),%xmm9
   0x0000000000002780 <+10112>:	addps  0x1d0(%r9),%xmm8
   0x0000000000002794 <+10132>:	addps  0x1e0(%r9),%xmm12
   0x00000000000027cf <+10191>:	addps  0x1f0(%r9),%xmm0

900	        _mm_store_ps(&C[i*4+C_OFFSET_24], C_row[i]);
   0x0000000000002774 <+10100>:	movaps %xmm9,0x1c0(%r9)
   0x0000000000002788 <+10120>:	movaps %xmm8,0x1d0(%r9)
   0x000000000000279c <+10140>:	movaps %xmm12,0x1e0(%r9)
   0x00000000000027d7 <+10199>:	movaps %xmm0,0x1f0(%r9)

901	      }
902	    }
903
904	    /* Reset C(3,1) matrix accumulators */
905	    C_row[0] = _mm_setzero_ps();
906	    C_row[1] = _mm_setzero_ps();
907	    C_row[2] = _mm_setzero_ps();
908	    C_row[3] = _mm_setzero_ps();
909
910	    if (norm[8]*norm[16] >= tolerance &&
   0x00000000000027df <+10207>:	movss  0x38(%rdx,%rsi,1),%xmm0
   0x00000000000027e5 <+10213>:	movaps %xmm6,%xmm7
   0x00000000000027e8 <+10216>:	mulss  %xmm0,%xmm7
   0x00000000000027ec <+10220>:	comiss %xmm1,%xmm7
   0x00000000000027ef <+10223>:	jb     0x2ce8 <stream_kernel+11496>

911	        norm[9]*norm[20] >= tolerance &&
   0x00000000000027f5 <+10229>:	movss  0x3c(%rdx,%rsi,1),%xmm7
   0x00000000000027fb <+10235>:	mulss  0x68(%rdx,%rsi,1),%xmm7
   0x0000000000002801 <+10241>:	comiss %xmm1,%xmm7
   0x0000000000002804 <+10244>:	jb     0x2ce8 <stream_kernel+11496>

912	        norm[10]*norm[24] >= tolerance &&
   0x000000000000280a <+10250>:	movss  0x40(%rdx,%rsi,1),%xmm7
   0x0000000000002810 <+10256>:	mulss  0x78(%rdx,%rsi,1),%xmm7
   0x0000000000002816 <+10262>:	comiss %xmm1,%xmm7
   0x0000000000002819 <+10265>:	jb     0x2ce8 <stream_kernel+11496>

913	        norm[11]*norm[28] >= tolerance)
   0x000000000000281f <+10271>:	movss  0x44(%rdx,%rsi,1),%xmm7
   0x0000000000002825 <+10277>:	mulss  0x88(%rdx,%rsi,1),%xmm7
   0x000000000000282e <+10286>:	comiss %xmm1,%xmm7
   0x0000000000002831 <+10289>:	jb     0x2ce8 <stream_kernel+11496>

914	    {
915	      /* A(3,1)*B(1,1) = C(3,1). */
916	      for (i = 0; i < 4; i++)
917	      {
918	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
919	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
   0x0000000000002837 <+10295>:	movaps (%rcx),%xmm7

920	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002849 <+10313>:	movaps 0x800(%r8),%xmm9
   0x0000000000002869 <+10345>:	movaps 0x840(%r8),%xmm15
   0x0000000000002871 <+10353>:	mulps  %xmm7,%xmm9
   0x000000000000287d <+10365>:	movaps 0x880(%r8),%xmm9
   0x00000000000028a5 <+10405>:	mulps  %xmm7,%xmm15
   0x00000000000028d9 <+10457>:	mulps  %xmm7,%xmm9
   0x00000000000028e5 <+10469>:	movaps 0x8c0(%r8),%xmm9
   0x000000000000290d <+10509>:	mulps  %xmm7,%xmm9

921
922	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
923	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
924	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000283a <+10298>:	movaps 0x10(%rcx),%xmm11
   0x0000000000002851 <+10321>:	movaps 0x810(%r8),%xmm13
   0x0000000000002875 <+10357>:	mulps  %xmm11,%xmm13
   0x0000000000002879 <+10361>:	addps  %xmm9,%xmm13
   0x000000000000288d <+10381>:	movaps 0x850(%r8),%xmm13
   0x00000000000028a9 <+10409>:	mulps  %xmm11,%xmm13
   0x00000000000028ad <+10413>:	addps  %xmm15,%xmm13
   0x00000000000028d1 <+10449>:	movaps 0x890(%r8),%xmm14
   0x00000000000028dd <+10461>:	mulps  %xmm11,%xmm14
   0x00000000000028e1 <+10465>:	addps  %xmm9,%xmm14
   0x0000000000002911 <+10513>:	movaps 0x8d0(%r8),%xmm7
   0x0000000000002919 <+10521>:	mulps  %xmm11,%xmm7
   0x0000000000002925 <+10533>:	addps  %xmm9,%xmm7

925
926	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
927	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
928	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000283f <+10303>:	movaps 0x20(%rcx),%xmm10
   0x0000000000002859 <+10329>:	movaps 0x820(%r8),%xmm14
   0x0000000000002885 <+10373>:	mulps  %xmm10,%xmm14
   0x0000000000002889 <+10377>:	addps  %xmm13,%xmm14
   0x000000000000289d <+10397>:	movaps 0x860(%r8),%xmm14
   0x00000000000028b1 <+10417>:	movaps 0x8a0(%r8),%xmm15
   0x00000000000028b9 <+10425>:	mulps  %xmm10,%xmm14
   0x00000000000028bd <+10429>:	addps  %xmm13,%xmm14
   0x00000000000028ed <+10477>:	mulps  %xmm10,%xmm15
   0x00000000000028f1 <+10481>:	addps  %xmm14,%xmm15
   0x000000000000291d <+10525>:	movaps 0x8e0(%r8),%xmm11
   0x0000000000002929 <+10537>:	mulps  %xmm10,%xmm11
   0x0000000000002935 <+10549>:	addps  %xmm7,%xmm11

929
930	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
931	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
932	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002844 <+10308>:	movaps 0x30(%rcx),%xmm12
   0x0000000000002861 <+10337>:	movaps 0x830(%r8),%xmm8
   0x0000000000002895 <+10389>:	mulps  %xmm12,%xmm8
   0x0000000000002899 <+10393>:	addps  %xmm14,%xmm8
   0x00000000000028c1 <+10433>:	movaps 0x870(%r8),%xmm13
   0x00000000000028c9 <+10441>:	mulps  %xmm12,%xmm13
   0x00000000000028cd <+10445>:	addps  %xmm14,%xmm13
   0x00000000000028f5 <+10485>:	movaps 0x8b0(%r8),%xmm14
   0x00000000000028fd <+10493>:	mulps  %xmm12,%xmm14
   0x0000000000002901 <+10497>:	addps  %xmm15,%xmm14
   0x0000000000002939 <+10553>:	movaps 0x8f0(%r8),%xmm7
   0x0000000000002941 <+10561>:	mulps  %xmm12,%xmm7
   0x000000000000294d <+10573>:	addps  %xmm11,%xmm7

933	      }
934
935	      /* A(3,2)*B(2,1) = C(3,1). */
936	      for (i = 0; i < 4; i++)
937	      {
938	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
939	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
940	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002945 <+10565>:	movaps 0x900(%r8),%xmm12
   0x0000000000002951 <+10577>:	movaps 0x100(%rcx),%xmm9
   0x0000000000002961 <+10593>:	mulps  %xmm9,%xmm12
   0x0000000000002965 <+10597>:	addps  %xmm8,%xmm12
   0x0000000000002999 <+10649>:	movaps 0x940(%r8),%xmm15
   0x00000000000029a1 <+10657>:	mulps  %xmm9,%xmm15
   0x00000000000029a5 <+10661>:	addps  %xmm13,%xmm15
   0x00000000000029d9 <+10713>:	movaps 0x980(%r8),%xmm15
   0x00000000000029e1 <+10721>:	mulps  %xmm9,%xmm15
   0x00000000000029e5 <+10725>:	addps  %xmm14,%xmm15
   0x0000000000002a19 <+10777>:	movaps 0x9c0(%r8),%xmm15
   0x0000000000002a21 <+10785>:	mulps  %xmm9,%xmm15
   0x0000000000002a31 <+10801>:	addps  %xmm7,%xmm15

941
942	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
943	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
944	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000292d <+10541>:	movaps 0x910(%r8),%xmm10
   0x0000000000002959 <+10585>:	movaps 0x110(%rcx),%xmm11
   0x0000000000002971 <+10609>:	mulps  %xmm11,%xmm10
   0x0000000000002975 <+10613>:	addps  %xmm12,%xmm10
   0x00000000000029a9 <+10665>:	movaps 0x950(%r8),%xmm13
   0x00000000000029b1 <+10673>:	mulps  %xmm11,%xmm13
   0x00000000000029b5 <+10677>:	addps  %xmm15,%xmm13
   0x00000000000029e9 <+10729>:	movaps 0x990(%r8),%xmm14
   0x00000000000029f1 <+10737>:	mulps  %xmm11,%xmm14
   0x00000000000029f5 <+10741>:	addps  %xmm15,%xmm14
   0x0000000000002a25 <+10789>:	movaps 0x9d0(%r8),%xmm9
   0x0000000000002a2d <+10797>:	mulps  %xmm11,%xmm9
   0x0000000000002a49 <+10825>:	addps  %xmm15,%xmm9

945
946	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
947	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
948	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002905 <+10501>:	movaps 0x920(%r8),%xmm15
   0x0000000000002969 <+10601>:	movaps 0x120(%rcx),%xmm8
   0x0000000000002981 <+10625>:	mulps  %xmm8,%xmm15
   0x0000000000002985 <+10629>:	addps  %xmm10,%xmm15
   0x00000000000029b9 <+10681>:	movaps 0x960(%r8),%xmm15
   0x00000000000029c1 <+10689>:	mulps  %xmm8,%xmm15
   0x00000000000029c5 <+10693>:	addps  %xmm13,%xmm15
   0x00000000000029f9 <+10745>:	movaps 0x9a0(%r8),%xmm15
   0x0000000000002a01 <+10753>:	mulps  %xmm8,%xmm15
   0x0000000000002a05 <+10757>:	addps  %xmm14,%xmm15
   0x0000000000002a4d <+10829>:	movaps 0x9e0(%r8),%xmm15
   0x0000000000002a55 <+10837>:	mulps  %xmm8,%xmm15
   0x0000000000002a59 <+10841>:	addps  %xmm9,%xmm15

949
950	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
951	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
952	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002979 <+10617>:	movaps 0x130(%rcx),%xmm12
   0x0000000000002989 <+10633>:	movaps 0x930(%r8),%xmm10
   0x0000000000002991 <+10641>:	mulps  %xmm12,%xmm10
   0x0000000000002995 <+10645>:	addps  %xmm15,%xmm10
   0x00000000000029c9 <+10697>:	movaps 0x970(%r8),%xmm13
   0x00000000000029d1 <+10705>:	mulps  %xmm12,%xmm13
   0x00000000000029d5 <+10709>:	addps  %xmm15,%xmm13
   0x0000000000002a09 <+10761>:	movaps 0x9b0(%r8),%xmm14
   0x0000000000002a11 <+10769>:	mulps  %xmm12,%xmm14
   0x0000000000002a15 <+10773>:	addps  %xmm15,%xmm14
   0x0000000000002a35 <+10805>:	movaps 0x9f0(%r8),%xmm7
   0x0000000000002a3d <+10813>:	mulps  %xmm12,%xmm7
   0x0000000000002a79 <+10873>:	addps  %xmm15,%xmm7

953	      }
954
955	      /* A(3,3)*B(3,1) = C(3,1). */
956	      for (i = 0; i < 4; i++)
957	      {
958	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
959	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
960	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002a41 <+10817>:	movaps 0xa00(%r8),%xmm12
   0x0000000000002a5d <+10845>:	movaps 0x200(%rcx),%xmm11
   0x0000000000002a75 <+10869>:	mulps  %xmm11,%xmm12
   0x0000000000002a89 <+10889>:	addps  %xmm10,%xmm12
   0x0000000000002ab9 <+10937>:	movaps 0xa40(%r8),%xmm15
   0x0000000000002ac1 <+10945>:	mulps  %xmm11,%xmm15
   0x0000000000002ac5 <+10949>:	addps  %xmm13,%xmm15
   0x0000000000002af9 <+11001>:	movaps 0xa80(%r8),%xmm15
   0x0000000000002b01 <+11009>:	mulps  %xmm11,%xmm15
   0x0000000000002b05 <+11013>:	addps  %xmm14,%xmm15
   0x0000000000002b39 <+11065>:	movaps 0xac0(%r8),%xmm15
   0x0000000000002b41 <+11073>:	mulps  %xmm11,%xmm15
   0x0000000000002b51 <+11089>:	addps  %xmm7,%xmm15

961
962	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
963	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
964	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002a65 <+10853>:	movaps 0x210(%rcx),%xmm8
   0x0000000000002a8d <+10893>:	movaps 0xa10(%r8),%xmm10
   0x0000000000002a95 <+10901>:	mulps  %xmm8,%xmm10
   0x0000000000002a99 <+10905>:	addps  %xmm12,%xmm10
   0x0000000000002ac9 <+10953>:	movaps 0xa50(%r8),%xmm13
   0x0000000000002ad1 <+10961>:	mulps  %xmm8,%xmm13
   0x0000000000002ad5 <+10965>:	addps  %xmm15,%xmm13
   0x0000000000002b09 <+11017>:	movaps 0xa90(%r8),%xmm14
   0x0000000000002b11 <+11025>:	mulps  %xmm8,%xmm14
   0x0000000000002b15 <+11029>:	addps  %xmm15,%xmm14
   0x0000000000002b55 <+11093>:	movaps 0xad0(%r8),%xmm7
   0x0000000000002b5d <+11101>:	mulps  %xmm8,%xmm7
   0x0000000000002b75 <+11125>:	addps  %xmm15,%xmm7

965
966	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
967	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
968	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002a6d <+10861>:	movaps 0x220(%rcx),%xmm9
   0x0000000000002a7d <+10877>:	movaps 0xa20(%r8),%xmm15
   0x0000000000002a85 <+10885>:	mulps  %xmm9,%xmm15
   0x0000000000002aa5 <+10917>:	addps  %xmm10,%xmm15
   0x0000000000002ad9 <+10969>:	movaps 0xa60(%r8),%xmm15
   0x0000000000002ae1 <+10977>:	mulps  %xmm9,%xmm15
   0x0000000000002ae5 <+10981>:	addps  %xmm13,%xmm15
   0x0000000000002b19 <+11033>:	movaps 0xaa0(%r8),%xmm15
   0x0000000000002b21 <+11041>:	mulps  %xmm9,%xmm15
   0x0000000000002b25 <+11045>:	addps  %xmm14,%xmm15
   0x0000000000002b45 <+11077>:	movaps 0xae0(%r8),%xmm11
   0x0000000000002b4d <+11085>:	mulps  %xmm9,%xmm11
   0x0000000000002b81 <+11137>:	addps  %xmm7,%xmm11

969
970	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
971	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
972	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002a9d <+10909>:	movaps 0xa30(%r8),%xmm12
   0x0000000000002aa9 <+10921>:	movaps 0x230(%rcx),%xmm10
   0x0000000000002ab1 <+10929>:	mulps  %xmm10,%xmm12
   0x0000000000002ab5 <+10933>:	addps  %xmm15,%xmm12
   0x0000000000002ae9 <+10985>:	movaps 0xa70(%r8),%xmm13
   0x0000000000002af1 <+10993>:	mulps  %xmm10,%xmm13
   0x0000000000002af5 <+10997>:	addps  %xmm15,%xmm13
   0x0000000000002b29 <+11049>:	movaps 0xab0(%r8),%xmm14
   0x0000000000002b31 <+11057>:	mulps  %xmm10,%xmm14
   0x0000000000002b35 <+11061>:	addps  %xmm15,%xmm14
   0x0000000000002b61 <+11105>:	movaps 0xaf0(%r8),%xmm8
   0x0000000000002b69 <+11113>:	mulps  %xmm10,%xmm8
   0x0000000000002b98 <+11160>:	addps  %xmm11,%xmm8

973	      }
974
975	      /* A(3,4)*B(4,1) = C(3,1). */
976	      for (i = 0; i < 4; i++)
977	      {
978	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
979	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
980	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b6d <+11117>:	movaps 0xb00(%r8),%xmm10
   0x0000000000002b85 <+11141>:	movaps 0x300(%rcx),%xmm9
   0x0000000000002b94 <+11156>:	mulps  %xmm9,%xmm10
   0x0000000000002ba8 <+11176>:	addps  %xmm12,%xmm10
   0x0000000000002bdc <+11228>:	movaps 0xb40(%r8),%xmm15
   0x0000000000002bf4 <+11252>:	mulps  %xmm9,%xmm15
   0x0000000000002bf8 <+11256>:	addps  %xmm13,%xmm15
   0x0000000000002c2c <+11308>:	movaps 0xb80(%r8),%xmm15
   0x0000000000002c34 <+11316>:	mulps  %xmm9,%xmm15
   0x0000000000002c38 <+11320>:	addps  %xmm14,%xmm15
   0x0000000000002c6c <+11372>:	movaps 0xbc0(%r8),%xmm15
   0x0000000000002c74 <+11380>:	mulps  %xmm9,%xmm15
   0x0000000000002c80 <+11392>:	addps  %xmm8,%xmm15

981
982	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
983	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
984	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b8d <+11149>:	movaps 0x310(%rcx),%xmm7
   0x0000000000002bac <+11180>:	movaps 0xb10(%r8),%xmm12
   0x0000000000002bb4 <+11188>:	mulps  %xmm7,%xmm12
   0x0000000000002bb8 <+11192>:	addps  %xmm10,%xmm12
   0x0000000000002bfc <+11260>:	movaps 0xb50(%r8),%xmm13
   0x0000000000002c04 <+11268>:	mulps  %xmm7,%xmm13
   0x0000000000002c08 <+11272>:	addps  %xmm15,%xmm13
   0x0000000000002c3c <+11324>:	movaps 0xb90(%r8),%xmm14
   0x0000000000002c44 <+11332>:	mulps  %xmm7,%xmm14
   0x0000000000002c48 <+11336>:	addps  %xmm15,%xmm14
   0x0000000000002c84 <+11396>:	movaps 0xbd0(%r8),%xmm8
   0x0000000000002c8c <+11404>:	mulps  %xmm7,%xmm8
   0x0000000000002c9c <+11420>:	addps  %xmm15,%xmm8

985
986	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
987	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
988	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b79 <+11129>:	movaps 0xb20(%r8),%xmm15
   0x0000000000002b9c <+11164>:	movaps 0x320(%rcx),%xmm11
   0x0000000000002ba4 <+11172>:	mulps  %xmm11,%xmm15
   0x0000000000002bc4 <+11204>:	addps  %xmm12,%xmm15
   0x0000000000002c0c <+11276>:	movaps 0xb60(%r8),%xmm15
   0x0000000000002c14 <+11284>:	mulps  %xmm11,%xmm15
   0x0000000000002c18 <+11288>:	addps  %xmm13,%xmm15
   0x0000000000002c4c <+11340>:	movaps 0xba0(%r8),%xmm15
   0x0000000000002c54 <+11348>:	mulps  %xmm11,%xmm15
   0x0000000000002c58 <+11352>:	addps  %xmm14,%xmm15
   0x0000000000002c90 <+11408>:	movaps 0xbe0(%r8),%xmm7
   0x0000000000002c98 <+11416>:	mulps  %xmm11,%xmm7
   0x0000000000002ca4 <+11428>:	addps  %xmm8,%xmm7

989
990	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
991	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
992	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002bbc <+11196>:	movaps 0xb30(%r8),%xmm10
   0x0000000000002bc8 <+11208>:	movaps 0x330(%rcx),%xmm12
   0x0000000000002bd0 <+11216>:	mulps  %xmm12,%xmm10
   0x0000000000002bd4 <+11220>:	addps  %xmm15,%xmm10
   0x0000000000002c1c <+11292>:	movaps 0xb70(%r8),%xmm13
   0x0000000000002c24 <+11300>:	mulps  %xmm12,%xmm13
   0x0000000000002c28 <+11304>:	addps  %xmm15,%xmm13
   0x0000000000002c5c <+11356>:	movaps 0xbb0(%r8),%xmm14
   0x0000000000002c64 <+11364>:	mulps  %xmm12,%xmm14
   0x0000000000002c68 <+11368>:	addps  %xmm15,%xmm14
   0x0000000000002c78 <+11384>:	movaps 0xbf0(%r8),%xmm9
   0x0000000000002ca0 <+11424>:	mulps  %xmm12,%xmm9
   0x0000000000002cd0 <+11472>:	addps  %xmm7,%xmm9

993	      }
994
995	      /* Store C(3,1) block. */
996	      for (i = 0; i < 4; i++)
997	      {
998	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000002bd8 <+11224>:	mulps  %xmm2,%xmm10
   0x0000000000002ca8 <+11432>:	mulps  %xmm2,%xmm13
   0x0000000000002cbc <+11452>:	mulps  %xmm2,%xmm14
   0x0000000000002cd4 <+11476>:	mulps  %xmm2,%xmm9

999	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_31]), C_row[i]);
   0x0000000000002be4 <+11236>:	addps  0x200(%r9),%xmm10
   0x0000000000002cac <+11436>:	addps  0x210(%r9),%xmm13
   0x0000000000002cc0 <+11456>:	addps  0x220(%r9),%xmm14
   0x0000000000002cd8 <+11480>:	addps  0x230(%r9),%xmm9

1000	        _mm_store_ps(&C[i*4+C_OFFSET_31], C_row[i]);
   0x0000000000002bec <+11244>:	movaps %xmm10,0x200(%r9)
   0x0000000000002cb4 <+11444>:	movaps %xmm13,0x210(%r9)
   0x0000000000002cc8 <+11464>:	movaps %xmm14,0x220(%r9)
   0x0000000000002ce0 <+11488>:	movaps %xmm9,0x230(%r9)

1001	      }
1002	    }
1003
1004	    /* Reset C(3,2) matrix accumulators */
1005	    C_row[0] = _mm_setzero_ps();
1006	    C_row[1] = _mm_setzero_ps();
1007	    C_row[2] = _mm_setzero_ps();
1008	    C_row[3] = _mm_setzero_ps();
1009
1010	    if (norm[8]*norm[17] >= tolerance &&
   0x0000000000002ce8 <+11496>:	movaps %xmm5,%xmm7
   0x0000000000002ceb <+11499>:	mulss  %xmm0,%xmm7
   0x0000000000002cef <+11503>:	comiss %xmm1,%xmm7
   0x0000000000002cf2 <+11506>:	jb     0x31ec <stream_kernel+12780>

1011	        norm[9]*norm[21] >= tolerance &&
   0x0000000000002cf8 <+11512>:	movss  0x3c(%rdx,%rsi,1),%xmm7
   0x0000000000002cfe <+11518>:	mulss  0x6c(%rdx,%rsi,1),%xmm7
   0x0000000000002d04 <+11524>:	comiss %xmm1,%xmm7
   0x0000000000002d07 <+11527>:	jb     0x31ec <stream_kernel+12780>

1012	        norm[10]*norm[25] >= tolerance &&
   0x0000000000002d0d <+11533>:	movss  0x40(%rdx,%rsi,1),%xmm7
   0x0000000000002d13 <+11539>:	mulss  0x7c(%rdx,%rsi,1),%xmm7
   0x0000000000002d19 <+11545>:	comiss %xmm1,%xmm7
   0x0000000000002d1c <+11548>:	jb     0x31ec <stream_kernel+12780>

1013	        norm[11]*norm[29] >= tolerance)
   0x0000000000002d22 <+11554>:	movss  0x44(%rdx,%rsi,1),%xmm7
   0x0000000000002d28 <+11560>:	mulss  0x8c(%rdx,%rsi,1),%xmm7
   0x0000000000002d31 <+11569>:	comiss %xmm1,%xmm7
   0x0000000000002d34 <+11572>:	jb     0x31ec <stream_kernel+12780>

1014	    {
1015	      /* A(3,1)*B(1,2) = C(3,2). */
1016	      for (i = 0; i < 4; i++)
1017	      {
1018	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
1019	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
1020	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002d3a <+11578>:	movaps 0x40(%rcx),%xmm7
   0x0000000000002d4d <+11597>:	movaps 0x800(%r8),%xmm9
   0x0000000000002d6d <+11629>:	movaps 0x840(%r8),%xmm15
   0x0000000000002d75 <+11637>:	mulps  %xmm7,%xmm9
   0x0000000000002d81 <+11649>:	movaps 0x880(%r8),%xmm9
   0x0000000000002da9 <+11689>:	mulps  %xmm7,%xmm15
   0x0000000000002ddd <+11741>:	mulps  %xmm7,%xmm9
   0x0000000000002de9 <+11753>:	movaps 0x8c0(%r8),%xmm9
   0x0000000000002e11 <+11793>:	mulps  %xmm7,%xmm9

1021
1022	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
1023	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
1024	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002d3e <+11582>:	movaps 0x50(%rcx),%xmm11
   0x0000000000002d55 <+11605>:	movaps 0x810(%r8),%xmm13
   0x0000000000002d79 <+11641>:	mulps  %xmm11,%xmm13
   0x0000000000002d7d <+11645>:	addps  %xmm9,%xmm13
   0x0000000000002d91 <+11665>:	movaps 0x850(%r8),%xmm13
   0x0000000000002dad <+11693>:	mulps  %xmm11,%xmm13
   0x0000000000002db1 <+11697>:	addps  %xmm15,%xmm13
   0x0000000000002dd5 <+11733>:	movaps 0x890(%r8),%xmm14
   0x0000000000002de1 <+11745>:	mulps  %xmm11,%xmm14
   0x0000000000002de5 <+11749>:	addps  %xmm9,%xmm14
   0x0000000000002e15 <+11797>:	movaps 0x8d0(%r8),%xmm7
   0x0000000000002e1d <+11805>:	mulps  %xmm11,%xmm7
   0x0000000000002e29 <+11817>:	addps  %xmm9,%xmm7

1025
1026	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
1027	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
1028	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002d43 <+11587>:	movaps 0x60(%rcx),%xmm10
   0x0000000000002d5d <+11613>:	movaps 0x820(%r8),%xmm14
   0x0000000000002d89 <+11657>:	mulps  %xmm10,%xmm14
   0x0000000000002d8d <+11661>:	addps  %xmm13,%xmm14
   0x0000000000002da1 <+11681>:	movaps 0x860(%r8),%xmm14
   0x0000000000002db5 <+11701>:	movaps 0x8a0(%r8),%xmm15
   0x0000000000002dbd <+11709>:	mulps  %xmm10,%xmm14
   0x0000000000002dc1 <+11713>:	addps  %xmm13,%xmm14
   0x0000000000002df1 <+11761>:	mulps  %xmm10,%xmm15
   0x0000000000002df5 <+11765>:	addps  %xmm14,%xmm15
   0x0000000000002e21 <+11809>:	movaps 0x8e0(%r8),%xmm11
   0x0000000000002e2d <+11821>:	mulps  %xmm10,%xmm11
   0x0000000000002e39 <+11833>:	addps  %xmm7,%xmm11

1029
1030	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
1031	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
1032	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002d48 <+11592>:	movaps 0x70(%rcx),%xmm12
   0x0000000000002d65 <+11621>:	movaps 0x830(%r8),%xmm8
   0x0000000000002d99 <+11673>:	mulps  %xmm12,%xmm8
   0x0000000000002d9d <+11677>:	addps  %xmm14,%xmm8
   0x0000000000002dc5 <+11717>:	movaps 0x870(%r8),%xmm13
   0x0000000000002dcd <+11725>:	mulps  %xmm12,%xmm13
   0x0000000000002dd1 <+11729>:	addps  %xmm14,%xmm13
   0x0000000000002df9 <+11769>:	movaps 0x8b0(%r8),%xmm14
   0x0000000000002e01 <+11777>:	mulps  %xmm12,%xmm14
   0x0000000000002e05 <+11781>:	addps  %xmm15,%xmm14
   0x0000000000002e3d <+11837>:	movaps 0x8f0(%r8),%xmm7
   0x0000000000002e45 <+11845>:	mulps  %xmm12,%xmm7
   0x0000000000002e51 <+11857>:	addps  %xmm11,%xmm7

1033	      }
1034
1035	      /* A(3,2)*B(2,2) = C(3,2). */
1036	      for (i = 0; i < 4; i++)
1037	      {
1038	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1039	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
1040	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002e49 <+11849>:	movaps 0x900(%r8),%xmm12
   0x0000000000002e55 <+11861>:	movaps 0x140(%rcx),%xmm9
   0x0000000000002e65 <+11877>:	mulps  %xmm9,%xmm12
   0x0000000000002e69 <+11881>:	addps  %xmm8,%xmm12
   0x0000000000002e9d <+11933>:	movaps 0x940(%r8),%xmm15
   0x0000000000002ea5 <+11941>:	mulps  %xmm9,%xmm15
   0x0000000000002ea9 <+11945>:	addps  %xmm13,%xmm15
   0x0000000000002edd <+11997>:	movaps 0x980(%r8),%xmm15
   0x0000000000002ee5 <+12005>:	mulps  %xmm9,%xmm15
   0x0000000000002ee9 <+12009>:	addps  %xmm14,%xmm15
   0x0000000000002f1d <+12061>:	movaps 0x9c0(%r8),%xmm15
   0x0000000000002f25 <+12069>:	mulps  %xmm9,%xmm15
   0x0000000000002f35 <+12085>:	addps  %xmm7,%xmm15

1041
1042	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1043	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
1044	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002e31 <+11825>:	movaps 0x910(%r8),%xmm10
   0x0000000000002e5d <+11869>:	movaps 0x150(%rcx),%xmm11
   0x0000000000002e75 <+11893>:	mulps  %xmm11,%xmm10
   0x0000000000002e79 <+11897>:	addps  %xmm12,%xmm10
   0x0000000000002ead <+11949>:	movaps 0x950(%r8),%xmm13
   0x0000000000002eb5 <+11957>:	mulps  %xmm11,%xmm13
   0x0000000000002eb9 <+11961>:	addps  %xmm15,%xmm13
   0x0000000000002eed <+12013>:	movaps 0x990(%r8),%xmm14
   0x0000000000002ef5 <+12021>:	mulps  %xmm11,%xmm14
   0x0000000000002ef9 <+12025>:	addps  %xmm15,%xmm14
   0x0000000000002f29 <+12073>:	movaps 0x9d0(%r8),%xmm9
   0x0000000000002f31 <+12081>:	mulps  %xmm11,%xmm9
   0x0000000000002f4d <+12109>:	addps  %xmm15,%xmm9

1045
1046	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1047	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
1048	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002e09 <+11785>:	movaps 0x920(%r8),%xmm15
   0x0000000000002e6d <+11885>:	movaps 0x160(%rcx),%xmm8
   0x0000000000002e85 <+11909>:	mulps  %xmm8,%xmm15
   0x0000000000002e89 <+11913>:	addps  %xmm10,%xmm15
   0x0000000000002ebd <+11965>:	movaps 0x960(%r8),%xmm15
   0x0000000000002ec5 <+11973>:	mulps  %xmm8,%xmm15
   0x0000000000002ec9 <+11977>:	addps  %xmm13,%xmm15
   0x0000000000002efd <+12029>:	movaps 0x9a0(%r8),%xmm15
   0x0000000000002f05 <+12037>:	mulps  %xmm8,%xmm15
   0x0000000000002f09 <+12041>:	addps  %xmm14,%xmm15
   0x0000000000002f51 <+12113>:	movaps 0x9e0(%r8),%xmm15
   0x0000000000002f59 <+12121>:	mulps  %xmm8,%xmm15
   0x0000000000002f5d <+12125>:	addps  %xmm9,%xmm15

1049
1050	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1051	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
1052	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002e7d <+11901>:	movaps 0x170(%rcx),%xmm12
   0x0000000000002e8d <+11917>:	movaps 0x930(%r8),%xmm10
   0x0000000000002e95 <+11925>:	mulps  %xmm12,%xmm10
   0x0000000000002e99 <+11929>:	addps  %xmm15,%xmm10
   0x0000000000002ecd <+11981>:	movaps 0x970(%r8),%xmm13
   0x0000000000002ed5 <+11989>:	mulps  %xmm12,%xmm13
   0x0000000000002ed9 <+11993>:	addps  %xmm15,%xmm13
   0x0000000000002f0d <+12045>:	movaps 0x9b0(%r8),%xmm14
   0x0000000000002f15 <+12053>:	mulps  %xmm12,%xmm14
   0x0000000000002f19 <+12057>:	addps  %xmm15,%xmm14
   0x0000000000002f39 <+12089>:	movaps 0x9f0(%r8),%xmm7
   0x0000000000002f41 <+12097>:	mulps  %xmm12,%xmm7
   0x0000000000002f7d <+12157>:	addps  %xmm15,%xmm7

1053	      }
1054
1055	      /* A(3,3)*B(3,2) = C(3,2). */
1056	      for (i = 0; i < 4; i++)
1057	      {
1058	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1059	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
1060	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002f45 <+12101>:	movaps 0xa00(%r8),%xmm12
   0x0000000000002f61 <+12129>:	movaps 0x240(%rcx),%xmm11
   0x0000000000002f79 <+12153>:	mulps  %xmm11,%xmm12
   0x0000000000002f8d <+12173>:	addps  %xmm10,%xmm12
   0x0000000000002fbd <+12221>:	movaps 0xa40(%r8),%xmm15
   0x0000000000002fc5 <+12229>:	mulps  %xmm11,%xmm15
   0x0000000000002fc9 <+12233>:	addps  %xmm13,%xmm15
   0x0000000000002ffd <+12285>:	movaps 0xa80(%r8),%xmm15
   0x0000000000003005 <+12293>:	mulps  %xmm11,%xmm15
   0x0000000000003009 <+12297>:	addps  %xmm14,%xmm15
   0x000000000000303d <+12349>:	movaps 0xac0(%r8),%xmm15
   0x0000000000003045 <+12357>:	mulps  %xmm11,%xmm15
   0x0000000000003055 <+12373>:	addps  %xmm7,%xmm15

1061
1062	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1063	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
1064	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002f69 <+12137>:	movaps 0x250(%rcx),%xmm8
   0x0000000000002f91 <+12177>:	movaps 0xa10(%r8),%xmm10
   0x0000000000002f99 <+12185>:	mulps  %xmm8,%xmm10
   0x0000000000002f9d <+12189>:	addps  %xmm12,%xmm10
   0x0000000000002fcd <+12237>:	movaps 0xa50(%r8),%xmm13
   0x0000000000002fd5 <+12245>:	mulps  %xmm8,%xmm13
   0x0000000000002fd9 <+12249>:	addps  %xmm15,%xmm13
   0x000000000000300d <+12301>:	movaps 0xa90(%r8),%xmm14
   0x0000000000003015 <+12309>:	mulps  %xmm8,%xmm14
   0x0000000000003019 <+12313>:	addps  %xmm15,%xmm14
   0x0000000000003059 <+12377>:	movaps 0xad0(%r8),%xmm7
   0x0000000000003061 <+12385>:	mulps  %xmm8,%xmm7
   0x0000000000003079 <+12409>:	addps  %xmm15,%xmm7

1065
1066	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1067	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
1068	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002f71 <+12145>:	movaps 0x260(%rcx),%xmm9
   0x0000000000002f81 <+12161>:	movaps 0xa20(%r8),%xmm15
   0x0000000000002f89 <+12169>:	mulps  %xmm9,%xmm15
   0x0000000000002fa9 <+12201>:	addps  %xmm10,%xmm15
   0x0000000000002fdd <+12253>:	movaps 0xa60(%r8),%xmm15
   0x0000000000002fe5 <+12261>:	mulps  %xmm9,%xmm15
   0x0000000000002fe9 <+12265>:	addps  %xmm13,%xmm15
   0x000000000000301d <+12317>:	movaps 0xaa0(%r8),%xmm15
   0x0000000000003025 <+12325>:	mulps  %xmm9,%xmm15
   0x0000000000003029 <+12329>:	addps  %xmm14,%xmm15
   0x0000000000003049 <+12361>:	movaps 0xae0(%r8),%xmm11
   0x0000000000003051 <+12369>:	mulps  %xmm9,%xmm11
   0x0000000000003085 <+12421>:	addps  %xmm7,%xmm11

1069
1070	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1071	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
1072	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002fa1 <+12193>:	movaps 0xa30(%r8),%xmm12
   0x0000000000002fad <+12205>:	movaps 0x270(%rcx),%xmm10
   0x0000000000002fb5 <+12213>:	mulps  %xmm10,%xmm12
   0x0000000000002fb9 <+12217>:	addps  %xmm15,%xmm12
   0x0000000000002fed <+12269>:	movaps 0xa70(%r8),%xmm13
   0x0000000000002ff5 <+12277>:	mulps  %xmm10,%xmm13
   0x0000000000002ff9 <+12281>:	addps  %xmm15,%xmm13
   0x000000000000302d <+12333>:	movaps 0xab0(%r8),%xmm14
   0x0000000000003035 <+12341>:	mulps  %xmm10,%xmm14
   0x0000000000003039 <+12345>:	addps  %xmm15,%xmm14
   0x0000000000003065 <+12389>:	movaps 0xaf0(%r8),%xmm8
   0x000000000000306d <+12397>:	mulps  %xmm10,%xmm8
   0x000000000000309c <+12444>:	addps  %xmm11,%xmm8

1073	      }
1074
1075	      /* A(3,4)*B(4,2) = C(3,2). */
1076	      for (i = 0; i < 4; i++)
1077	      {
1078	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1079	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
1080	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003071 <+12401>:	movaps 0xb00(%r8),%xmm10
   0x0000000000003089 <+12425>:	movaps 0x340(%rcx),%xmm9
   0x0000000000003098 <+12440>:	mulps  %xmm9,%xmm10
   0x00000000000030ac <+12460>:	addps  %xmm12,%xmm10
   0x00000000000030e0 <+12512>:	movaps 0xb40(%r8),%xmm15
   0x00000000000030f8 <+12536>:	mulps  %xmm9,%xmm15
   0x00000000000030fc <+12540>:	addps  %xmm13,%xmm15
   0x0000000000003130 <+12592>:	movaps 0xb80(%r8),%xmm15
   0x0000000000003138 <+12600>:	mulps  %xmm9,%xmm15
   0x000000000000313c <+12604>:	addps  %xmm14,%xmm15
   0x0000000000003170 <+12656>:	movaps 0xbc0(%r8),%xmm15
   0x0000000000003178 <+12664>:	mulps  %xmm9,%xmm15
   0x0000000000003184 <+12676>:	addps  %xmm8,%xmm15

1081
1082	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1083	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
1084	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003091 <+12433>:	movaps 0x350(%rcx),%xmm7
   0x00000000000030b0 <+12464>:	movaps 0xb10(%r8),%xmm12
   0x00000000000030b8 <+12472>:	mulps  %xmm7,%xmm12
   0x00000000000030bc <+12476>:	addps  %xmm10,%xmm12
   0x0000000000003100 <+12544>:	movaps 0xb50(%r8),%xmm13
   0x0000000000003108 <+12552>:	mulps  %xmm7,%xmm13
   0x000000000000310c <+12556>:	addps  %xmm15,%xmm13
   0x0000000000003140 <+12608>:	movaps 0xb90(%r8),%xmm14
   0x0000000000003148 <+12616>:	mulps  %xmm7,%xmm14
   0x000000000000314c <+12620>:	addps  %xmm15,%xmm14
   0x0000000000003188 <+12680>:	movaps 0xbd0(%r8),%xmm8
   0x0000000000003190 <+12688>:	mulps  %xmm7,%xmm8
   0x00000000000031a0 <+12704>:	addps  %xmm15,%xmm8

1085
1086	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1087	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
1088	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000307d <+12413>:	movaps 0xb20(%r8),%xmm15
   0x00000000000030a0 <+12448>:	movaps 0x360(%rcx),%xmm11
   0x00000000000030a8 <+12456>:	mulps  %xmm11,%xmm15
   0x00000000000030c8 <+12488>:	addps  %xmm12,%xmm15
   0x0000000000003110 <+12560>:	movaps 0xb60(%r8),%xmm15
   0x0000000000003118 <+12568>:	mulps  %xmm11,%xmm15
   0x000000000000311c <+12572>:	addps  %xmm13,%xmm15
   0x0000000000003150 <+12624>:	movaps 0xba0(%r8),%xmm15
   0x0000000000003158 <+12632>:	mulps  %xmm11,%xmm15
   0x000000000000315c <+12636>:	addps  %xmm14,%xmm15
   0x0000000000003194 <+12692>:	movaps 0xbe0(%r8),%xmm7
   0x000000000000319c <+12700>:	mulps  %xmm11,%xmm7
   0x00000000000031a8 <+12712>:	addps  %xmm8,%xmm7

1089
1090	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1091	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
1092	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000030c0 <+12480>:	movaps 0xb30(%r8),%xmm10
   0x00000000000030cc <+12492>:	movaps 0x370(%rcx),%xmm12
   0x00000000000030d4 <+12500>:	mulps  %xmm12,%xmm10
   0x00000000000030d8 <+12504>:	addps  %xmm15,%xmm10
   0x0000000000003120 <+12576>:	movaps 0xb70(%r8),%xmm13
   0x0000000000003128 <+12584>:	mulps  %xmm12,%xmm13
   0x000000000000312c <+12588>:	addps  %xmm15,%xmm13
   0x0000000000003160 <+12640>:	movaps 0xbb0(%r8),%xmm14
   0x0000000000003168 <+12648>:	mulps  %xmm12,%xmm14
   0x000000000000316c <+12652>:	addps  %xmm15,%xmm14
   0x000000000000317c <+12668>:	movaps 0xbf0(%r8),%xmm9
   0x00000000000031a4 <+12708>:	mulps  %xmm12,%xmm9
   0x00000000000031d4 <+12756>:	addps  %xmm7,%xmm9

1093	      }
1094
1095	      /* Store C(3,2) block. */
1096	      for (i = 0; i < 4; i++)
1097	      {
1098	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000030dc <+12508>:	mulps  %xmm2,%xmm10
   0x00000000000031ac <+12716>:	mulps  %xmm2,%xmm13
   0x00000000000031c0 <+12736>:	mulps  %xmm2,%xmm14
   0x00000000000031d8 <+12760>:	mulps  %xmm2,%xmm9

1099	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_32]), C_row[i]);
   0x00000000000030e8 <+12520>:	addps  0x240(%r9),%xmm10
   0x00000000000031b0 <+12720>:	addps  0x250(%r9),%xmm13
   0x00000000000031c4 <+12740>:	addps  0x260(%r9),%xmm14
   0x00000000000031dc <+12764>:	addps  0x270(%r9),%xmm9

1100	        _mm_store_ps(&C[i*4+C_OFFSET_32], C_row[i]);
   0x00000000000030f0 <+12528>:	movaps %xmm10,0x240(%r9)
   0x00000000000031b8 <+12728>:	movaps %xmm13,0x250(%r9)
   0x00000000000031cc <+12748>:	movaps %xmm14,0x260(%r9)
   0x00000000000031e4 <+12772>:	movaps %xmm9,0x270(%r9)

1101	      }
1102	    }
1103
1104	    /* Reset C(3,3) matrix accumulators */
1105	    C_row[0] = _mm_setzero_ps();
1106	    C_row[1] = _mm_setzero_ps();
1107	    C_row[2] = _mm_setzero_ps();
1108	    C_row[3] = _mm_setzero_ps();
1109
1110	    if (norm[8]*norm[18] >= tolerance &&
   0x00000000000031ec <+12780>:	movaps %xmm4,%xmm7
   0x00000000000031ef <+12783>:	mulss  %xmm0,%xmm7
   0x00000000000031f3 <+12787>:	comiss %xmm1,%xmm7
   0x00000000000031f6 <+12790>:	jb     0x36ff <stream_kernel+14079>

1111	        norm[9]*norm[22] >= tolerance &&
   0x00000000000031fc <+12796>:	movss  0x3c(%rdx,%rsi,1),%xmm7
   0x0000000000003202 <+12802>:	mulss  0x70(%rdx,%rsi,1),%xmm7
   0x0000000000003208 <+12808>:	comiss %xmm1,%xmm7
   0x000000000000320b <+12811>:	jb     0x36ff <stream_kernel+14079>

1112	        norm[10]*norm[26] >= tolerance &&
   0x0000000000003211 <+12817>:	movss  0x40(%rdx,%rsi,1),%xmm7
   0x0000000000003217 <+12823>:	mulss  0x80(%rdx,%rsi,1),%xmm7
   0x0000000000003220 <+12832>:	comiss %xmm1,%xmm7
   0x0000000000003223 <+12835>:	jb     0x36ff <stream_kernel+14079>

1113	        norm[11]*norm[30] >= tolerance)
   0x0000000000003229 <+12841>:	movss  0x44(%rdx,%rsi,1),%xmm7
   0x000000000000322f <+12847>:	mulss  0x90(%rdx,%rsi,1),%xmm7
   0x0000000000003238 <+12856>:	comiss %xmm1,%xmm7
   0x000000000000323b <+12859>:	jb     0x36ff <stream_kernel+14079>

1114	    {
1115	      /* A(3,1)*B(1,3) = C(3,3). */
1116	      for (i = 0; i < 4; i++)
1117	      {
1118	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
1119	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
1120	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003241 <+12865>:	movaps 0x80(%rcx),%xmm7
   0x0000000000003260 <+12896>:	movaps 0x800(%r8),%xmm9
   0x0000000000003280 <+12928>:	movaps 0x840(%r8),%xmm15
   0x0000000000003288 <+12936>:	mulps  %xmm7,%xmm9
   0x0000000000003294 <+12948>:	movaps 0x880(%r8),%xmm9
   0x00000000000032bc <+12988>:	mulps  %xmm7,%xmm15
   0x00000000000032f0 <+13040>:	mulps  %xmm7,%xmm9
   0x00000000000032fc <+13052>:	movaps 0x8c0(%r8),%xmm9
   0x0000000000003324 <+13092>:	mulps  %xmm7,%xmm9

1121
1122	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
1123	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
1124	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003248 <+12872>:	movaps 0x90(%rcx),%xmm11
   0x0000000000003268 <+12904>:	movaps 0x810(%r8),%xmm13
   0x000000000000328c <+12940>:	mulps  %xmm11,%xmm13
   0x0000000000003290 <+12944>:	addps  %xmm9,%xmm13
   0x00000000000032a4 <+12964>:	movaps 0x850(%r8),%xmm13
   0x00000000000032c0 <+12992>:	mulps  %xmm11,%xmm13
   0x00000000000032c4 <+12996>:	addps  %xmm15,%xmm13
   0x00000000000032e8 <+13032>:	movaps 0x890(%r8),%xmm14
   0x00000000000032f4 <+13044>:	mulps  %xmm11,%xmm14
   0x00000000000032f8 <+13048>:	addps  %xmm9,%xmm14
   0x0000000000003328 <+13096>:	movaps 0x8d0(%r8),%xmm7
   0x0000000000003330 <+13104>:	mulps  %xmm11,%xmm7
   0x000000000000333c <+13116>:	addps  %xmm9,%xmm7

1125
1126	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
1127	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
1128	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003250 <+12880>:	movaps 0xa0(%rcx),%xmm10
   0x0000000000003270 <+12912>:	movaps 0x820(%r8),%xmm14
   0x000000000000329c <+12956>:	mulps  %xmm10,%xmm14
   0x00000000000032a0 <+12960>:	addps  %xmm13,%xmm14
   0x00000000000032b4 <+12980>:	movaps 0x860(%r8),%xmm14
   0x00000000000032c8 <+13000>:	movaps 0x8a0(%r8),%xmm15
   0x00000000000032d0 <+13008>:	mulps  %xmm10,%xmm14
   0x00000000000032d4 <+13012>:	addps  %xmm13,%xmm14
   0x0000000000003304 <+13060>:	mulps  %xmm10,%xmm15
   0x0000000000003308 <+13064>:	addps  %xmm14,%xmm15
   0x0000000000003334 <+13108>:	movaps 0x8e0(%r8),%xmm11
   0x0000000000003340 <+13120>:	mulps  %xmm10,%xmm11
   0x000000000000334c <+13132>:	addps  %xmm7,%xmm11

1129
1130	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
1131	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
1132	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003258 <+12888>:	movaps 0xb0(%rcx),%xmm12
   0x0000000000003278 <+12920>:	movaps 0x830(%r8),%xmm8
   0x00000000000032ac <+12972>:	mulps  %xmm12,%xmm8
   0x00000000000032b0 <+12976>:	addps  %xmm14,%xmm8
   0x00000000000032d8 <+13016>:	movaps 0x870(%r8),%xmm13
   0x00000000000032e0 <+13024>:	mulps  %xmm12,%xmm13
   0x00000000000032e4 <+13028>:	addps  %xmm14,%xmm13
   0x000000000000330c <+13068>:	movaps 0x8b0(%r8),%xmm14
   0x0000000000003314 <+13076>:	mulps  %xmm12,%xmm14
   0x0000000000003318 <+13080>:	addps  %xmm15,%xmm14
   0x0000000000003350 <+13136>:	movaps 0x8f0(%r8),%xmm7
   0x0000000000003358 <+13144>:	mulps  %xmm12,%xmm7
   0x0000000000003364 <+13156>:	addps  %xmm11,%xmm7

1133	      }
1134
1135	      /* A(3,2)*B(2,3) = C(3,3). */
1136	      for (i = 0; i < 4; i++)
1137	      {
1138	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1139	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
1140	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000335c <+13148>:	movaps 0x900(%r8),%xmm12
   0x0000000000003368 <+13160>:	movaps 0x180(%rcx),%xmm9
   0x0000000000003378 <+13176>:	mulps  %xmm9,%xmm12
   0x000000000000337c <+13180>:	addps  %xmm8,%xmm12
   0x00000000000033b0 <+13232>:	movaps 0x940(%r8),%xmm15
   0x00000000000033b8 <+13240>:	mulps  %xmm9,%xmm15
   0x00000000000033bc <+13244>:	addps  %xmm13,%xmm15
   0x00000000000033f0 <+13296>:	movaps 0x980(%r8),%xmm15
   0x00000000000033f8 <+13304>:	mulps  %xmm9,%xmm15
   0x00000000000033fc <+13308>:	addps  %xmm14,%xmm15
   0x0000000000003430 <+13360>:	movaps 0x9c0(%r8),%xmm15
   0x0000000000003438 <+13368>:	mulps  %xmm9,%xmm15
   0x0000000000003448 <+13384>:	addps  %xmm7,%xmm15

1141
1142	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1143	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
1144	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003344 <+13124>:	movaps 0x910(%r8),%xmm10
   0x0000000000003370 <+13168>:	movaps 0x190(%rcx),%xmm11
   0x0000000000003388 <+13192>:	mulps  %xmm11,%xmm10
   0x000000000000338c <+13196>:	addps  %xmm12,%xmm10
   0x00000000000033c0 <+13248>:	movaps 0x950(%r8),%xmm13
   0x00000000000033c8 <+13256>:	mulps  %xmm11,%xmm13
   0x00000000000033cc <+13260>:	addps  %xmm15,%xmm13
   0x0000000000003400 <+13312>:	movaps 0x990(%r8),%xmm14
   0x0000000000003408 <+13320>:	mulps  %xmm11,%xmm14
   0x000000000000340c <+13324>:	addps  %xmm15,%xmm14
   0x000000000000343c <+13372>:	movaps 0x9d0(%r8),%xmm9
   0x0000000000003444 <+13380>:	mulps  %xmm11,%xmm9
   0x0000000000003460 <+13408>:	addps  %xmm15,%xmm9

1145
1146	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1147	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
1148	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000331c <+13084>:	movaps 0x920(%r8),%xmm15
   0x0000000000003380 <+13184>:	movaps 0x1a0(%rcx),%xmm8
   0x0000000000003398 <+13208>:	mulps  %xmm8,%xmm15
   0x000000000000339c <+13212>:	addps  %xmm10,%xmm15
   0x00000000000033d0 <+13264>:	movaps 0x960(%r8),%xmm15
   0x00000000000033d8 <+13272>:	mulps  %xmm8,%xmm15
   0x00000000000033dc <+13276>:	addps  %xmm13,%xmm15
   0x0000000000003410 <+13328>:	movaps 0x9a0(%r8),%xmm15
   0x0000000000003418 <+13336>:	mulps  %xmm8,%xmm15
   0x000000000000341c <+13340>:	addps  %xmm14,%xmm15
   0x0000000000003464 <+13412>:	movaps 0x9e0(%r8),%xmm15
   0x000000000000346c <+13420>:	mulps  %xmm8,%xmm15
   0x0000000000003470 <+13424>:	addps  %xmm9,%xmm15

1149
1150	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1151	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
1152	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003390 <+13200>:	movaps 0x1b0(%rcx),%xmm12
   0x00000000000033a0 <+13216>:	movaps 0x930(%r8),%xmm10
   0x00000000000033a8 <+13224>:	mulps  %xmm12,%xmm10
   0x00000000000033ac <+13228>:	addps  %xmm15,%xmm10
   0x00000000000033e0 <+13280>:	movaps 0x970(%r8),%xmm13
   0x00000000000033e8 <+13288>:	mulps  %xmm12,%xmm13
   0x00000000000033ec <+13292>:	addps  %xmm15,%xmm13
   0x0000000000003420 <+13344>:	movaps 0x9b0(%r8),%xmm14
   0x0000000000003428 <+13352>:	mulps  %xmm12,%xmm14
   0x000000000000342c <+13356>:	addps  %xmm15,%xmm14
   0x000000000000344c <+13388>:	movaps 0x9f0(%r8),%xmm7
   0x0000000000003454 <+13396>:	mulps  %xmm12,%xmm7
   0x0000000000003490 <+13456>:	addps  %xmm15,%xmm7

1153	      }
1154
1155	      /* A(3,3)*B(3,3) = C(3,3). */
1156	      for (i = 0; i < 4; i++)
1157	      {
1158	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1159	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
1160	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003458 <+13400>:	movaps 0xa00(%r8),%xmm12
   0x0000000000003474 <+13428>:	movaps 0x280(%rcx),%xmm11
   0x000000000000348c <+13452>:	mulps  %xmm11,%xmm12
   0x00000000000034a0 <+13472>:	addps  %xmm10,%xmm12
   0x00000000000034d0 <+13520>:	movaps 0xa40(%r8),%xmm15
   0x00000000000034d8 <+13528>:	mulps  %xmm11,%xmm15
   0x00000000000034dc <+13532>:	addps  %xmm13,%xmm15
   0x0000000000003510 <+13584>:	movaps 0xa80(%r8),%xmm15
   0x0000000000003518 <+13592>:	mulps  %xmm11,%xmm15
   0x000000000000351c <+13596>:	addps  %xmm14,%xmm15
   0x0000000000003550 <+13648>:	movaps 0xac0(%r8),%xmm15
   0x0000000000003558 <+13656>:	mulps  %xmm11,%xmm15
   0x0000000000003568 <+13672>:	addps  %xmm7,%xmm15

1161
1162	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1163	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
1164	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000347c <+13436>:	movaps 0x290(%rcx),%xmm8
   0x00000000000034a4 <+13476>:	movaps 0xa10(%r8),%xmm10
   0x00000000000034ac <+13484>:	mulps  %xmm8,%xmm10
   0x00000000000034b0 <+13488>:	addps  %xmm12,%xmm10
   0x00000000000034e0 <+13536>:	movaps 0xa50(%r8),%xmm13
   0x00000000000034e8 <+13544>:	mulps  %xmm8,%xmm13
   0x00000000000034ec <+13548>:	addps  %xmm15,%xmm13
   0x0000000000003520 <+13600>:	movaps 0xa90(%r8),%xmm14
   0x0000000000003528 <+13608>:	mulps  %xmm8,%xmm14
   0x000000000000352c <+13612>:	addps  %xmm15,%xmm14
   0x000000000000356c <+13676>:	movaps 0xad0(%r8),%xmm7
   0x0000000000003574 <+13684>:	mulps  %xmm8,%xmm7
   0x000000000000358c <+13708>:	addps  %xmm15,%xmm7

1165
1166	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1167	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
1168	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003484 <+13444>:	movaps 0x2a0(%rcx),%xmm9
   0x0000000000003494 <+13460>:	movaps 0xa20(%r8),%xmm15
   0x000000000000349c <+13468>:	mulps  %xmm9,%xmm15
   0x00000000000034bc <+13500>:	addps  %xmm10,%xmm15
   0x00000000000034f0 <+13552>:	movaps 0xa60(%r8),%xmm15
   0x00000000000034f8 <+13560>:	mulps  %xmm9,%xmm15
   0x00000000000034fc <+13564>:	addps  %xmm13,%xmm15
   0x0000000000003530 <+13616>:	movaps 0xaa0(%r8),%xmm15
   0x0000000000003538 <+13624>:	mulps  %xmm9,%xmm15
   0x000000000000353c <+13628>:	addps  %xmm14,%xmm15
   0x000000000000355c <+13660>:	movaps 0xae0(%r8),%xmm11
   0x0000000000003564 <+13668>:	mulps  %xmm9,%xmm11
   0x0000000000003598 <+13720>:	addps  %xmm7,%xmm11

1169
1170	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1171	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
1172	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000034b4 <+13492>:	movaps 0xa30(%r8),%xmm12
   0x00000000000034c0 <+13504>:	movaps 0x2b0(%rcx),%xmm10
   0x00000000000034c8 <+13512>:	mulps  %xmm10,%xmm12
   0x00000000000034cc <+13516>:	addps  %xmm15,%xmm12
   0x0000000000003500 <+13568>:	movaps 0xa70(%r8),%xmm13
   0x0000000000003508 <+13576>:	mulps  %xmm10,%xmm13
   0x000000000000350c <+13580>:	addps  %xmm15,%xmm13
   0x0000000000003540 <+13632>:	movaps 0xab0(%r8),%xmm14
   0x0000000000003548 <+13640>:	mulps  %xmm10,%xmm14
   0x000000000000354c <+13644>:	addps  %xmm15,%xmm14
   0x0000000000003578 <+13688>:	movaps 0xaf0(%r8),%xmm8
   0x0000000000003580 <+13696>:	mulps  %xmm10,%xmm8
   0x00000000000035af <+13743>:	addps  %xmm11,%xmm8

1173	      }
1174
1175	      /* A(3,4)*B(4,3) = C(3,3). */
1176	      for (i = 0; i < 4; i++)
1177	      {
1178	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1179	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
1180	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003584 <+13700>:	movaps 0xb00(%r8),%xmm10
   0x000000000000359c <+13724>:	movaps 0x380(%rcx),%xmm9
   0x00000000000035ab <+13739>:	mulps  %xmm9,%xmm10
   0x00000000000035bf <+13759>:	addps  %xmm12,%xmm10
   0x00000000000035f3 <+13811>:	movaps 0xb40(%r8),%xmm15
   0x000000000000360b <+13835>:	mulps  %xmm9,%xmm15
   0x000000000000360f <+13839>:	addps  %xmm13,%xmm15
   0x0000000000003643 <+13891>:	movaps 0xb80(%r8),%xmm15
   0x000000000000364b <+13899>:	mulps  %xmm9,%xmm15
   0x000000000000364f <+13903>:	addps  %xmm14,%xmm15
   0x0000000000003683 <+13955>:	movaps 0xbc0(%r8),%xmm15
   0x000000000000368b <+13963>:	mulps  %xmm9,%xmm15
   0x0000000000003697 <+13975>:	addps  %xmm8,%xmm15

1181
1182	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1183	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
1184	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000035a4 <+13732>:	movaps 0x390(%rcx),%xmm7
   0x00000000000035c3 <+13763>:	movaps 0xb10(%r8),%xmm12
   0x00000000000035cb <+13771>:	mulps  %xmm7,%xmm12
   0x00000000000035cf <+13775>:	addps  %xmm10,%xmm12
   0x0000000000003613 <+13843>:	movaps 0xb50(%r8),%xmm13
   0x000000000000361b <+13851>:	mulps  %xmm7,%xmm13
   0x000000000000361f <+13855>:	addps  %xmm15,%xmm13
   0x0000000000003653 <+13907>:	movaps 0xb90(%r8),%xmm14
   0x000000000000365b <+13915>:	mulps  %xmm7,%xmm14
   0x000000000000365f <+13919>:	addps  %xmm15,%xmm14
   0x000000000000369b <+13979>:	movaps 0xbd0(%r8),%xmm8
   0x00000000000036a3 <+13987>:	mulps  %xmm7,%xmm8
   0x00000000000036b3 <+14003>:	addps  %xmm15,%xmm8

1185
1186	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1187	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
1188	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003590 <+13712>:	movaps 0xb20(%r8),%xmm15
   0x00000000000035b3 <+13747>:	movaps 0x3a0(%rcx),%xmm11
   0x00000000000035bb <+13755>:	mulps  %xmm11,%xmm15
   0x00000000000035db <+13787>:	addps  %xmm12,%xmm15
   0x0000000000003623 <+13859>:	movaps 0xb60(%r8),%xmm15
   0x000000000000362b <+13867>:	mulps  %xmm11,%xmm15
   0x000000000000362f <+13871>:	addps  %xmm13,%xmm15
   0x0000000000003663 <+13923>:	movaps 0xba0(%r8),%xmm15
   0x000000000000366b <+13931>:	mulps  %xmm11,%xmm15
   0x000000000000366f <+13935>:	addps  %xmm14,%xmm15
   0x00000000000036a7 <+13991>:	movaps 0xbe0(%r8),%xmm7
   0x00000000000036af <+13999>:	mulps  %xmm11,%xmm7
   0x00000000000036bb <+14011>:	addps  %xmm8,%xmm7

1189
1190	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1191	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
1192	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000035d3 <+13779>:	movaps 0xb30(%r8),%xmm10
   0x00000000000035df <+13791>:	movaps 0x3b0(%rcx),%xmm12
   0x00000000000035e7 <+13799>:	mulps  %xmm12,%xmm10
   0x00000000000035eb <+13803>:	addps  %xmm15,%xmm10
   0x0000000000003633 <+13875>:	movaps 0xb70(%r8),%xmm13
   0x000000000000363b <+13883>:	mulps  %xmm12,%xmm13
   0x000000000000363f <+13887>:	addps  %xmm15,%xmm13
   0x0000000000003673 <+13939>:	movaps 0xbb0(%r8),%xmm14
   0x000000000000367b <+13947>:	mulps  %xmm12,%xmm14
   0x000000000000367f <+13951>:	addps  %xmm15,%xmm14
   0x000000000000368f <+13967>:	movaps 0xbf0(%r8),%xmm9
   0x00000000000036b7 <+14007>:	mulps  %xmm12,%xmm9
   0x00000000000036e7 <+14055>:	addps  %xmm7,%xmm9

1193	      }
1194
1195	      /* Store C(3,3) block. */
1196	      for (i = 0; i < 4; i++)
1197	      {
1198	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000035ef <+13807>:	mulps  %xmm2,%xmm10
   0x00000000000036bf <+14015>:	mulps  %xmm2,%xmm13
   0x00000000000036d3 <+14035>:	mulps  %xmm2,%xmm14
   0x00000000000036eb <+14059>:	mulps  %xmm2,%xmm9

1199	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_33]), C_row[i]);
   0x00000000000035fb <+13819>:	addps  0x280(%r9),%xmm10
   0x00000000000036c3 <+14019>:	addps  0x290(%r9),%xmm13
   0x00000000000036d7 <+14039>:	addps  0x2a0(%r9),%xmm14
   0x00000000000036ef <+14063>:	addps  0x2b0(%r9),%xmm9

1200	        _mm_store_ps(&C[i*4+C_OFFSET_33], C_row[i]);
   0x0000000000003603 <+13827>:	movaps %xmm10,0x280(%r9)
   0x00000000000036cb <+14027>:	movaps %xmm13,0x290(%r9)
   0x00000000000036df <+14047>:	movaps %xmm14,0x2a0(%r9)
   0x00000000000036f7 <+14071>:	movaps %xmm9,0x2b0(%r9)

1201	      }
1202	    }
1203
1204	    /* Reset C(3,4) matrix accumulators */
1205	    C_row[0] = _mm_setzero_ps();
1206	    C_row[1] = _mm_setzero_ps();
1207	    C_row[2] = _mm_setzero_ps();
1208	    C_row[3] = _mm_setzero_ps();
1209
1210	    if (norm[8]*norm[19] >= tolerance &&
   0x00000000000036ff <+14079>:	mulss  %xmm3,%xmm0
   0x0000000000003703 <+14083>:	comiss %xmm1,%xmm0
   0x0000000000003706 <+14086>:	jb     0x3c08 <stream_kernel+15368>

1211	        norm[9]*norm[23] >= tolerance &&
   0x000000000000370c <+14092>:	movss  0x3c(%rdx,%rsi,1),%xmm0
   0x0000000000003712 <+14098>:	mulss  0x74(%rdx,%rsi,1),%xmm0
   0x0000000000003718 <+14104>:	comiss %xmm1,%xmm0
   0x000000000000371b <+14107>:	jb     0x3c08 <stream_kernel+15368>

1212	        norm[10]*norm[27] >= tolerance &&
   0x0000000000003721 <+14113>:	movss  0x40(%rdx,%rsi,1),%xmm0
   0x0000000000003727 <+14119>:	mulss  0x84(%rdx,%rsi,1),%xmm0
   0x0000000000003730 <+14128>:	comiss %xmm1,%xmm0
   0x0000000000003733 <+14131>:	jb     0x3c08 <stream_kernel+15368>

1213	        norm[11]*norm[31] >= tolerance)
   0x0000000000003739 <+14137>:	movss  0x44(%rdx,%rsi,1),%xmm0
   0x000000000000373f <+14143>:	mulss  0x94(%rdx,%rsi,1),%xmm0
   0x0000000000003748 <+14152>:	comiss %xmm1,%xmm0
   0x000000000000374b <+14155>:	jb     0x3c08 <stream_kernel+15368>

1214	    {
1215	      /* A(3,1)*B(1,4) = C(3,4). */
1216	      for (i = 0; i < 4; i++)
1217	      {
1218	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
1219	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
1220	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003751 <+14161>:	movaps 0xc0(%rcx),%xmm7
   0x000000000000376f <+14191>:	movaps 0x800(%r8),%xmm10
   0x000000000000378f <+14223>:	movaps 0x840(%r8),%xmm15
   0x000000000000379f <+14239>:	mulps  %xmm7,%xmm10
   0x00000000000037cb <+14283>:	movaps 0x880(%r8),%xmm11
   0x00000000000037d3 <+14291>:	mulps  %xmm7,%xmm15
   0x00000000000037ef <+14319>:	movaps 0x8c0(%r8),%xmm14
   0x0000000000003807 <+14343>:	mulps  %xmm7,%xmm11
   0x000000000000383b <+14395>:	mulps  %xmm7,%xmm14

1221
1222	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
1223	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
1224	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003758 <+14168>:	movaps 0xd0(%rcx),%xmm0
   0x0000000000003777 <+14199>:	movaps 0x810(%r8),%xmm9
   0x0000000000003797 <+14231>:	movaps 0x850(%r8),%xmm14
   0x00000000000037a3 <+14243>:	mulps  %xmm0,%xmm9
   0x00000000000037a7 <+14247>:	addps  %xmm10,%xmm9
   0x00000000000037d7 <+14295>:	mulps  %xmm0,%xmm14
   0x00000000000037db <+14299>:	addps  %xmm15,%xmm14
   0x00000000000037df <+14303>:	movaps 0x890(%r8),%xmm15
   0x000000000000380b <+14347>:	mulps  %xmm0,%xmm15
   0x000000000000380f <+14351>:	addps  %xmm11,%xmm15
   0x000000000000383f <+14399>:	movaps 0x8d0(%r8),%xmm7
   0x0000000000003847 <+14407>:	mulps  %xmm0,%xmm7
   0x0000000000003852 <+14418>:	addps  %xmm14,%xmm7

1225
1226	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
1227	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
1228	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000375f <+14175>:	movaps 0xe0(%rcx),%xmm8
   0x000000000000377f <+14207>:	movaps 0x820(%r8),%xmm11
   0x00000000000037ab <+14251>:	movaps 0x860(%r8),%xmm10
   0x00000000000037b3 <+14259>:	mulps  %xmm8,%xmm11
   0x00000000000037b7 <+14263>:	addps  %xmm9,%xmm11
   0x00000000000037e7 <+14311>:	mulps  %xmm8,%xmm10
   0x00000000000037eb <+14315>:	addps  %xmm14,%xmm10
   0x0000000000003813 <+14355>:	movaps 0x8a0(%r8),%xmm11
   0x000000000000381b <+14363>:	mulps  %xmm8,%xmm11
   0x000000000000381f <+14367>:	addps  %xmm15,%xmm11
   0x000000000000384a <+14410>:	movaps 0x8e0(%r8),%xmm0
   0x0000000000003856 <+14422>:	mulps  %xmm8,%xmm0
   0x000000000000386e <+14446>:	addps  %xmm7,%xmm0

1229
1230	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
1231	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
1232	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003767 <+14183>:	movaps 0xf0(%rcx),%xmm12
   0x0000000000003787 <+14215>:	movaps 0x830(%r8),%xmm13
   0x00000000000037bb <+14267>:	movaps 0x870(%r8),%xmm9
   0x00000000000037c3 <+14275>:	mulps  %xmm12,%xmm13
   0x00000000000037c7 <+14279>:	addps  %xmm11,%xmm13
   0x00000000000037f7 <+14327>:	mulps  %xmm12,%xmm9
   0x00000000000037fb <+14331>:	addps  %xmm10,%xmm9
   0x00000000000037ff <+14335>:	movaps 0x8b0(%r8),%xmm10
   0x000000000000382b <+14379>:	mulps  %xmm12,%xmm10
   0x000000000000382f <+14383>:	addps  %xmm11,%xmm10
   0x0000000000003833 <+14387>:	movaps 0x8f0(%r8),%xmm11
   0x000000000000386a <+14442>:	mulps  %xmm12,%xmm11
   0x0000000000003885 <+14469>:	addps  %xmm0,%xmm11

1233	      }
1234
1235	      /* A(3,2)*B(2,4) = C(3,4). */
1236	      for (i = 0; i < 4; i++)
1237	      {
1238	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1239	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
1240	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000385a <+14426>:	movaps 0x1c0(%rcx),%xmm8
   0x0000000000003871 <+14449>:	movaps 0x900(%r8),%xmm12
   0x0000000000003881 <+14465>:	mulps  %xmm8,%xmm12
   0x0000000000003889 <+14473>:	movaps 0x940(%r8),%xmm0
   0x0000000000003895 <+14485>:	addps  %xmm13,%xmm12
   0x00000000000038a9 <+14505>:	mulps  %xmm8,%xmm0
   0x00000000000038c5 <+14533>:	addps  %xmm9,%xmm0
   0x00000000000038f9 <+14585>:	movaps 0x980(%r8),%xmm9
   0x0000000000003901 <+14593>:	mulps  %xmm8,%xmm9
   0x0000000000003914 <+14612>:	addps  %xmm10,%xmm9
   0x0000000000003938 <+14648>:	movaps 0x9c0(%r8),%xmm10
   0x0000000000003940 <+14656>:	mulps  %xmm8,%xmm10
   0x0000000000003964 <+14692>:	addps  %xmm11,%xmm10

1241
1242	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1243	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
1244	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003862 <+14434>:	movaps 0x1d0(%rcx),%xmm14
   0x0000000000003899 <+14489>:	movaps 0x910(%r8),%xmm13
   0x00000000000038a1 <+14497>:	mulps  %xmm14,%xmm13
   0x00000000000038a5 <+14501>:	addps  %xmm12,%xmm13
   0x00000000000038c9 <+14537>:	movaps 0x950(%r8),%xmm9
   0x00000000000038d1 <+14545>:	mulps  %xmm14,%xmm9
   0x00000000000038e5 <+14565>:	addps  %xmm0,%xmm9
   0x0000000000003918 <+14616>:	movaps 0x990(%r8),%xmm10
   0x0000000000003920 <+14624>:	mulps  %xmm14,%xmm10
   0x0000000000003924 <+14628>:	addps  %xmm9,%xmm10
   0x0000000000003944 <+14660>:	movaps 0x9d0(%r8),%xmm8
   0x000000000000394c <+14668>:	mulps  %xmm14,%xmm8
   0x0000000000003980 <+14720>:	addps  %xmm10,%xmm8

1245
1246	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1247	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
1248	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003823 <+14371>:	movaps 0x1e0(%rcx),%xmm15
   0x0000000000003879 <+14457>:	movaps 0x920(%r8),%xmm7
   0x0000000000003891 <+14481>:	mulps  %xmm15,%xmm7
   0x00000000000038b5 <+14517>:	addps  %xmm13,%xmm7
   0x00000000000038e9 <+14569>:	movaps 0x960(%r8),%xmm0
   0x00000000000038f1 <+14577>:	mulps  %xmm15,%xmm0
   0x00000000000038f5 <+14581>:	addps  %xmm9,%xmm0
   0x0000000000003928 <+14632>:	movaps 0x9a0(%r8),%xmm9
   0x0000000000003930 <+14640>:	mulps  %xmm15,%xmm9
   0x0000000000003934 <+14644>:	addps  %xmm10,%xmm9
   0x0000000000003984 <+14724>:	movaps 0x9e0(%r8),%xmm10
   0x000000000000398c <+14732>:	mulps  %xmm15,%xmm10
   0x00000000000039a0 <+14752>:	addps  %xmm8,%xmm10

1249
1250	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1251	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
1252	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000038ad <+14509>:	movaps 0x1f0(%rcx),%xmm12
   0x00000000000038b9 <+14521>:	movaps 0x930(%r8),%xmm13
   0x00000000000038c1 <+14529>:	mulps  %xmm12,%xmm13
   0x00000000000038d5 <+14549>:	addps  %xmm7,%xmm13
   0x00000000000038d9 <+14553>:	movaps 0x970(%r8),%xmm7
   0x00000000000038e1 <+14561>:	mulps  %xmm12,%xmm7
   0x0000000000003905 <+14597>:	addps  %xmm0,%xmm7
   0x0000000000003908 <+14600>:	movaps 0x9b0(%r8),%xmm0
   0x0000000000003910 <+14608>:	mulps  %xmm12,%xmm0
   0x0000000000003950 <+14672>:	addps  %xmm9,%xmm0
   0x0000000000003968 <+14696>:	movaps 0x9f0(%r8),%xmm11
   0x0000000000003970 <+14704>:	mulps  %xmm12,%xmm11
   0x00000000000039ac <+14764>:	addps  %xmm10,%xmm11

1253	      }
1254
1255	      /* A(3,3)*B(3,4) = C(3,4). */
1256	      for (i = 0; i < 4; i++)
1257	      {
1258	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1259	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
1260	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003954 <+14676>:	movaps 0x2c0(%rcx),%xmm14
   0x0000000000003974 <+14708>:	movaps 0xa00(%r8),%xmm12
   0x000000000000397c <+14716>:	mulps  %xmm14,%xmm12
   0x0000000000003990 <+14736>:	addps  %xmm13,%xmm12
   0x00000000000039c8 <+14792>:	movaps 0xa40(%r8),%xmm12
   0x00000000000039d0 <+14800>:	mulps  %xmm14,%xmm12
   0x00000000000039e4 <+14820>:	addps  %xmm7,%xmm12
   0x0000000000003a08 <+14856>:	movaps 0xa80(%r8),%xmm12
   0x0000000000003a10 <+14864>:	mulps  %xmm14,%xmm12
   0x0000000000003a24 <+14884>:	addps  %xmm0,%xmm12
   0x0000000000003a48 <+14920>:	movaps 0xac0(%r8),%xmm12
   0x0000000000003a50 <+14928>:	mulps  %xmm14,%xmm12
   0x0000000000003a70 <+14960>:	addps  %xmm11,%xmm12

1261
1262	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1263	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
1264	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000395c <+14684>:	movaps 0x2d0(%rcx),%xmm9
   0x0000000000003994 <+14740>:	movaps 0xa10(%r8),%xmm13
   0x000000000000399c <+14748>:	mulps  %xmm9,%xmm13
   0x00000000000039c4 <+14788>:	addps  %xmm12,%xmm13
   0x00000000000039e8 <+14824>:	movaps 0xa50(%r8),%xmm7
   0x00000000000039f0 <+14832>:	mulps  %xmm9,%xmm7
   0x0000000000003a04 <+14852>:	addps  %xmm12,%xmm7
   0x0000000000003a28 <+14888>:	movaps 0xa90(%r8),%xmm0
   0x0000000000003a30 <+14896>:	mulps  %xmm9,%xmm0
   0x0000000000003a44 <+14916>:	addps  %xmm12,%xmm0
   0x0000000000003a74 <+14964>:	movaps 0xad0(%r8),%xmm11
   0x0000000000003a7c <+14972>:	mulps  %xmm9,%xmm11
   0x0000000000003aa0 <+15008>:	addps  %xmm12,%xmm11

1265
1266	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1267	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
1268	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000039a4 <+14756>:	movaps 0xa20(%r8),%xmm8
   0x00000000000039b0 <+14768>:	movaps 0x2e0(%rcx),%xmm10
   0x00000000000039c0 <+14784>:	mulps  %xmm10,%xmm8
   0x00000000000039d4 <+14804>:	addps  %xmm13,%xmm8
   0x00000000000039f8 <+14840>:	movaps 0xa60(%r8),%xmm8
   0x0000000000003a00 <+14848>:	mulps  %xmm10,%xmm8
   0x0000000000003a14 <+14868>:	addps  %xmm7,%xmm8
   0x0000000000003a38 <+14904>:	movaps 0xaa0(%r8),%xmm8
   0x0000000000003a40 <+14912>:	mulps  %xmm10,%xmm8
   0x0000000000003a60 <+14944>:	addps  %xmm0,%xmm8
   0x0000000000003a80 <+14976>:	movaps 0xae0(%r8),%xmm9
   0x0000000000003a88 <+14984>:	mulps  %xmm10,%xmm9
   0x0000000000003aac <+15020>:	addps  %xmm11,%xmm9

1269
1270	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1271	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
1272	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000039b8 <+14776>:	movaps 0x2f0(%rcx),%xmm15
   0x00000000000039d8 <+14808>:	movaps 0xa30(%r8),%xmm13
   0x00000000000039e0 <+14816>:	mulps  %xmm15,%xmm13
   0x00000000000039f4 <+14836>:	addps  %xmm8,%xmm13
   0x0000000000003a18 <+14872>:	movaps 0xa70(%r8),%xmm7
   0x0000000000003a20 <+14880>:	mulps  %xmm15,%xmm7
   0x0000000000003a34 <+14900>:	addps  %xmm8,%xmm7
   0x0000000000003a54 <+14932>:	movaps 0xaf0(%r8),%xmm14
   0x0000000000003a5c <+14940>:	mulps  %xmm15,%xmm14
   0x0000000000003a64 <+14948>:	movaps 0xab0(%r8),%xmm0
   0x0000000000003a6c <+14956>:	mulps  %xmm15,%xmm0
   0x0000000000003a94 <+14996>:	addps  %xmm8,%xmm0
   0x0000000000003ac4 <+15044>:	addps  %xmm9,%xmm14

1273	      }
1274
1275	      /* A(3,4)*B(4,4) = C(3,4). */
1276	      for (i = 0; i < 4; i++)
1277	      {
1278	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1279	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
1280	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003a8c <+14988>:	movaps 0xb00(%r8),%xmm10
   0x0000000000003ab0 <+15024>:	movaps 0x3c0(%rcx),%xmm15
   0x0000000000003ac0 <+15040>:	mulps  %xmm15,%xmm10
   0x0000000000003ad4 <+15060>:	addps  %xmm13,%xmm10
   0x0000000000003af8 <+15096>:	movaps 0xb40(%r8),%xmm8
   0x0000000000003b00 <+15104>:	mulps  %xmm15,%xmm8
   0x0000000000003b14 <+15124>:	addps  %xmm7,%xmm8
   0x0000000000003b38 <+15160>:	movaps 0xb80(%r8),%xmm7
   0x0000000000003b40 <+15168>:	mulps  %xmm15,%xmm7
   0x0000000000003b54 <+15188>:	addps  %xmm0,%xmm7
   0x0000000000003b75 <+15221>:	movaps 0xbc0(%r8),%xmm0
   0x0000000000003b7d <+15229>:	mulps  %xmm15,%xmm0
   0x0000000000003bcd <+15309>:	addps  %xmm14,%xmm0

1281
1282	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1283	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
1284	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003a98 <+15000>:	movaps 0xb10(%r8),%xmm8
   0x0000000000003ab8 <+15032>:	movaps 0x3d0(%rcx),%xmm11
   0x0000000000003ac8 <+15048>:	mulps  %xmm11,%xmm8
   0x0000000000003ae4 <+15076>:	addps  %xmm10,%xmm8
   0x0000000000003b18 <+15128>:	movaps 0xb50(%r8),%xmm7
   0x0000000000003b20 <+15136>:	mulps  %xmm11,%xmm7
   0x0000000000003b24 <+15140>:	addps  %xmm8,%xmm7
   0x0000000000003b57 <+15191>:	movaps 0xb90(%r8),%xmm0
   0x0000000000003b5f <+15199>:	mulps  %xmm11,%xmm0
   0x0000000000003b63 <+15203>:	addps  %xmm7,%xmm0
   0x0000000000003bd1 <+15313>:	movaps 0xbd0(%r8),%xmm14
   0x0000000000003bd9 <+15321>:	mulps  %xmm11,%xmm14
   0x0000000000003bdd <+15325>:	addps  %xmm0,%xmm14

1285
1286	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1287	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
1288	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003aa4 <+15012>:	movaps 0xb20(%r8),%xmm12
   0x0000000000003ad8 <+15064>:	movaps 0x3e0(%rcx),%xmm13
   0x0000000000003ae0 <+15072>:	mulps  %xmm13,%xmm12
   0x0000000000003af4 <+15092>:	addps  %xmm8,%xmm12
   0x0000000000003b08 <+15112>:	movaps 0xb60(%r8),%xmm12
   0x0000000000003b10 <+15120>:	mulps  %xmm13,%xmm12
   0x0000000000003b34 <+15156>:	addps  %xmm7,%xmm12
   0x0000000000003b66 <+15206>:	movaps 0xba0(%r8),%xmm7
   0x0000000000003b6e <+15214>:	mulps  %xmm13,%xmm7
   0x0000000000003b72 <+15218>:	addps  %xmm0,%xmm7
   0x0000000000003b81 <+15233>:	movaps 0xbe0(%r8),%xmm15
   0x0000000000003b89 <+15241>:	mulps  %xmm13,%xmm15
   0x0000000000003bed <+15341>:	addps  %xmm14,%xmm15

1289
1290	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1291	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
1292	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003acc <+15052>:	movaps 0xb30(%r8),%xmm9
   0x0000000000003ae8 <+15080>:	movaps 0x3f0(%rcx),%xmm10
   0x0000000000003af0 <+15088>:	mulps  %xmm10,%xmm9
   0x0000000000003b04 <+15108>:	addps  %xmm12,%xmm9
   0x0000000000003b28 <+15144>:	movaps 0xb70(%r8),%xmm8
   0x0000000000003b30 <+15152>:	mulps  %xmm10,%xmm8
   0x0000000000003b44 <+15172>:	addps  %xmm12,%xmm8
   0x0000000000003b48 <+15176>:	movaps 0xbb0(%r8),%xmm12
   0x0000000000003b50 <+15184>:	mulps  %xmm10,%xmm12
   0x0000000000003b8d <+15245>:	addps  %xmm7,%xmm12
   0x0000000000003be1 <+15329>:	movaps 0xbf0(%r8),%xmm0
   0x0000000000003be9 <+15337>:	mulps  %xmm10,%xmm0
   0x0000000000003bf1 <+15345>:	addps  %xmm15,%xmm0

1293	      }
1294
1295	      /* Store C(3,4) block. */
1296	      for (i = 0; i < 4; i++)
1297	      {
1298	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000003b91 <+15249>:	mulps  %xmm2,%xmm9
   0x0000000000003ba5 <+15269>:	mulps  %xmm2,%xmm8
   0x0000000000003bb9 <+15289>:	mulps  %xmm2,%xmm12
   0x0000000000003bf5 <+15349>:	mulps  %xmm2,%xmm0

1299	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_34]), C_row[i]);
   0x0000000000003b95 <+15253>:	addps  0x2c0(%r9),%xmm9
   0x0000000000003ba9 <+15273>:	addps  0x2d0(%r9),%xmm8
   0x0000000000003bbd <+15293>:	addps  0x2e0(%r9),%xmm12
   0x0000000000003bf8 <+15352>:	addps  0x2f0(%r9),%xmm0

1300	        _mm_store_ps(&C[i*4+C_OFFSET_34], C_row[i]);
   0x0000000000003b9d <+15261>:	movaps %xmm9,0x2c0(%r9)
   0x0000000000003bb1 <+15281>:	movaps %xmm8,0x2d0(%r9)
   0x0000000000003bc5 <+15301>:	movaps %xmm12,0x2e0(%r9)
   0x0000000000003c00 <+15360>:	movaps %xmm0,0x2f0(%r9)

1301	      }
1302	    }
1303
1304	    /* Reset C(4,1) matrix accumulators */
1305	    C_row[0] = _mm_setzero_ps();
1306	    C_row[1] = _mm_setzero_ps();
1307	    C_row[2] = _mm_setzero_ps();
1308	    C_row[3] = _mm_setzero_ps();
1309
1310	    if (norm[12]*norm[16] >= tolerance &&
   0x0000000000003c08 <+15368>:	movss  0x48(%rdx,%rsi,1),%xmm7
   0x0000000000003c0e <+15374>:	mulss  %xmm7,%xmm6
   0x0000000000003c12 <+15378>:	comiss %xmm1,%xmm6
   0x0000000000003c15 <+15381>:	jb     0x4107 <stream_kernel+16647>

1311	        norm[13]*norm[20] >= tolerance &&
   0x0000000000003c1b <+15387>:	movss  0x4c(%rdx,%rsi,1),%xmm0
   0x0000000000003c21 <+15393>:	mulss  0x68(%rdx,%rsi,1),%xmm0
   0x0000000000003c27 <+15399>:	comiss %xmm1,%xmm0
   0x0000000000003c2a <+15402>:	jb     0x4107 <stream_kernel+16647>

1312	        norm[14]*norm[24] >= tolerance &&
   0x0000000000003c30 <+15408>:	movss  0x50(%rdx,%rsi,1),%xmm0
   0x0000000000003c36 <+15414>:	mulss  0x78(%rdx,%rsi,1),%xmm0
   0x0000000000003c3c <+15420>:	comiss %xmm1,%xmm0
   0x0000000000003c3f <+15423>:	jb     0x4107 <stream_kernel+16647>

1313	        norm[15]*norm[28] >= tolerance)
   0x0000000000003c45 <+15429>:	movss  0x54(%rdx,%rsi,1),%xmm0
   0x0000000000003c4b <+15435>:	mulss  0x88(%rdx,%rsi,1),%xmm0
   0x0000000000003c54 <+15444>:	comiss %xmm1,%xmm0
   0x0000000000003c57 <+15447>:	jb     0x4107 <stream_kernel+16647>

1314	    {
1315	      /* A(4,1)*B(1,1) = C(4,1). */
1316	      for (i = 0; i < 4; i++)
1317	      {
1318	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1319	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
   0x0000000000003c5d <+15453>:	movaps (%rcx),%xmm6

1320	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003c6e <+15470>:	movaps 0xc00(%r8),%xmm10
   0x0000000000003c8e <+15502>:	movaps 0xc40(%r8),%xmm15
   0x0000000000003c9e <+15518>:	mulps  %xmm6,%xmm10
   0x0000000000003cca <+15562>:	movaps 0xc80(%r8),%xmm11
   0x0000000000003cd2 <+15570>:	mulps  %xmm6,%xmm15
   0x0000000000003cee <+15598>:	movaps 0xcc0(%r8),%xmm14
   0x0000000000003d06 <+15622>:	mulps  %xmm6,%xmm11
   0x0000000000003d3a <+15674>:	mulps  %xmm6,%xmm14

1321
1322	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1323	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
1324	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003c60 <+15456>:	movaps 0x10(%rcx),%xmm0
   0x0000000000003c76 <+15478>:	movaps 0xc10(%r8),%xmm9
   0x0000000000003c96 <+15510>:	movaps 0xc50(%r8),%xmm14
   0x0000000000003ca2 <+15522>:	mulps  %xmm0,%xmm9
   0x0000000000003ca6 <+15526>:	addps  %xmm10,%xmm9
   0x0000000000003cd6 <+15574>:	mulps  %xmm0,%xmm14
   0x0000000000003cda <+15578>:	addps  %xmm15,%xmm14
   0x0000000000003cde <+15582>:	movaps 0xc90(%r8),%xmm15
   0x0000000000003d0a <+15626>:	mulps  %xmm0,%xmm15
   0x0000000000003d0e <+15630>:	addps  %xmm11,%xmm15
   0x0000000000003d3e <+15678>:	movaps 0xcd0(%r8),%xmm6
   0x0000000000003d46 <+15686>:	mulps  %xmm0,%xmm6
   0x0000000000003d51 <+15697>:	addps  %xmm14,%xmm6

1325
1326	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1327	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
1328	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003c64 <+15460>:	movaps 0x20(%rcx),%xmm8
   0x0000000000003c7e <+15486>:	movaps 0xc20(%r8),%xmm11
   0x0000000000003caa <+15530>:	movaps 0xc60(%r8),%xmm10
   0x0000000000003cb2 <+15538>:	mulps  %xmm8,%xmm11
   0x0000000000003cb6 <+15542>:	addps  %xmm9,%xmm11
   0x0000000000003ce6 <+15590>:	mulps  %xmm8,%xmm10
   0x0000000000003cea <+15594>:	addps  %xmm14,%xmm10
   0x0000000000003d12 <+15634>:	movaps 0xca0(%r8),%xmm11
   0x0000000000003d1a <+15642>:	mulps  %xmm8,%xmm11
   0x0000000000003d1e <+15646>:	addps  %xmm15,%xmm11
   0x0000000000003d49 <+15689>:	movaps 0xce0(%r8),%xmm0
   0x0000000000003d55 <+15701>:	mulps  %xmm8,%xmm0
   0x0000000000003d6d <+15725>:	addps  %xmm6,%xmm0

1329
1330	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1331	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
1332	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003c69 <+15465>:	movaps 0x30(%rcx),%xmm12
   0x0000000000003c86 <+15494>:	movaps 0xc30(%r8),%xmm13
   0x0000000000003cba <+15546>:	movaps 0xc70(%r8),%xmm9
   0x0000000000003cc2 <+15554>:	mulps  %xmm12,%xmm13
   0x0000000000003cc6 <+15558>:	addps  %xmm11,%xmm13
   0x0000000000003cf6 <+15606>:	mulps  %xmm12,%xmm9
   0x0000000000003cfa <+15610>:	addps  %xmm10,%xmm9
   0x0000000000003cfe <+15614>:	movaps 0xcb0(%r8),%xmm10
   0x0000000000003d2a <+15658>:	mulps  %xmm12,%xmm10
   0x0000000000003d2e <+15662>:	addps  %xmm11,%xmm10
   0x0000000000003d32 <+15666>:	movaps 0xcf0(%r8),%xmm11
   0x0000000000003d69 <+15721>:	mulps  %xmm12,%xmm11
   0x0000000000003d84 <+15748>:	addps  %xmm0,%xmm11

1333	      }
1334
1335	      /* A(4,2)*B(2,1) = C(4,1). */
1336	      for (i = 0; i < 4; i++)
1337	      {
1338	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1339	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
1340	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003d59 <+15705>:	movaps 0x100(%rcx),%xmm8
   0x0000000000003d70 <+15728>:	movaps 0xd00(%r8),%xmm12
   0x0000000000003d80 <+15744>:	mulps  %xmm8,%xmm12
   0x0000000000003d88 <+15752>:	movaps 0xd40(%r8),%xmm0
   0x0000000000003d94 <+15764>:	addps  %xmm13,%xmm12
   0x0000000000003da8 <+15784>:	mulps  %xmm8,%xmm0
   0x0000000000003dc4 <+15812>:	addps  %xmm9,%xmm0
   0x0000000000003df8 <+15864>:	movaps 0xd80(%r8),%xmm9
   0x0000000000003e00 <+15872>:	mulps  %xmm8,%xmm9
   0x0000000000003e13 <+15891>:	addps  %xmm10,%xmm9
   0x0000000000003e37 <+15927>:	movaps 0xdc0(%r8),%xmm10
   0x0000000000003e3f <+15935>:	mulps  %xmm8,%xmm10
   0x0000000000003e63 <+15971>:	addps  %xmm11,%xmm10

1341
1342	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1343	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
1344	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003d61 <+15713>:	movaps 0x110(%rcx),%xmm14
   0x0000000000003d98 <+15768>:	movaps 0xd10(%r8),%xmm13
   0x0000000000003da0 <+15776>:	mulps  %xmm14,%xmm13
   0x0000000000003da4 <+15780>:	addps  %xmm12,%xmm13
   0x0000000000003dc8 <+15816>:	movaps 0xd50(%r8),%xmm9
   0x0000000000003dd0 <+15824>:	mulps  %xmm14,%xmm9
   0x0000000000003de4 <+15844>:	addps  %xmm0,%xmm9
   0x0000000000003e17 <+15895>:	movaps 0xd90(%r8),%xmm10
   0x0000000000003e1f <+15903>:	mulps  %xmm14,%xmm10
   0x0000000000003e23 <+15907>:	addps  %xmm9,%xmm10
   0x0000000000003e43 <+15939>:	movaps 0xdd0(%r8),%xmm8
   0x0000000000003e4b <+15947>:	mulps  %xmm14,%xmm8
   0x0000000000003e7f <+15999>:	addps  %xmm10,%xmm8

1345
1346	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1347	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
1348	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003d22 <+15650>:	movaps 0x120(%rcx),%xmm15
   0x0000000000003d78 <+15736>:	movaps 0xd20(%r8),%xmm6
   0x0000000000003d90 <+15760>:	mulps  %xmm15,%xmm6
   0x0000000000003db4 <+15796>:	addps  %xmm13,%xmm6
   0x0000000000003de8 <+15848>:	movaps 0xd60(%r8),%xmm0
   0x0000000000003df0 <+15856>:	mulps  %xmm15,%xmm0
   0x0000000000003df4 <+15860>:	addps  %xmm9,%xmm0
   0x0000000000003e27 <+15911>:	movaps 0xda0(%r8),%xmm9
   0x0000000000003e2f <+15919>:	mulps  %xmm15,%xmm9
   0x0000000000003e33 <+15923>:	addps  %xmm10,%xmm9
   0x0000000000003e83 <+16003>:	movaps 0xde0(%r8),%xmm10
   0x0000000000003e8b <+16011>:	mulps  %xmm15,%xmm10
   0x0000000000003e9f <+16031>:	addps  %xmm8,%xmm10

1349
1350	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1351	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
1352	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003dac <+15788>:	movaps 0x130(%rcx),%xmm12
   0x0000000000003db8 <+15800>:	movaps 0xd30(%r8),%xmm13
   0x0000000000003dc0 <+15808>:	mulps  %xmm12,%xmm13
   0x0000000000003dd4 <+15828>:	addps  %xmm6,%xmm13
   0x0000000000003dd8 <+15832>:	movaps 0xd70(%r8),%xmm6
   0x0000000000003de0 <+15840>:	mulps  %xmm12,%xmm6
   0x0000000000003e04 <+15876>:	addps  %xmm0,%xmm6
   0x0000000000003e07 <+15879>:	movaps 0xdb0(%r8),%xmm0
   0x0000000000003e0f <+15887>:	mulps  %xmm12,%xmm0
   0x0000000000003e4f <+15951>:	addps  %xmm9,%xmm0
   0x0000000000003e67 <+15975>:	movaps 0xdf0(%r8),%xmm11
   0x0000000000003e6f <+15983>:	mulps  %xmm12,%xmm11
   0x0000000000003eab <+16043>:	addps  %xmm10,%xmm11

1353	      }
1354
1355	      /* A(4,3)*B(3,1) = C(4,1). */
1356	      for (i = 0; i < 4; i++)
1357	      {
1358	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1359	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
1360	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003e53 <+15955>:	movaps 0x200(%rcx),%xmm14
   0x0000000000003e73 <+15987>:	movaps 0xe00(%r8),%xmm12
   0x0000000000003e7b <+15995>:	mulps  %xmm14,%xmm12
   0x0000000000003e8f <+16015>:	addps  %xmm13,%xmm12
   0x0000000000003ec7 <+16071>:	movaps 0xe40(%r8),%xmm12
   0x0000000000003ecf <+16079>:	mulps  %xmm14,%xmm12
   0x0000000000003ee3 <+16099>:	addps  %xmm6,%xmm12
   0x0000000000003f07 <+16135>:	movaps 0xe80(%r8),%xmm12
   0x0000000000003f0f <+16143>:	mulps  %xmm14,%xmm12
   0x0000000000003f23 <+16163>:	addps  %xmm0,%xmm12
   0x0000000000003f47 <+16199>:	movaps 0xec0(%r8),%xmm12
   0x0000000000003f4f <+16207>:	mulps  %xmm14,%xmm12
   0x0000000000003f6f <+16239>:	addps  %xmm11,%xmm12

1361
1362	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1363	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
1364	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003e5b <+15963>:	movaps 0x210(%rcx),%xmm9
   0x0000000000003e93 <+16019>:	movaps 0xe10(%r8),%xmm13
   0x0000000000003e9b <+16027>:	mulps  %xmm9,%xmm13
   0x0000000000003ec3 <+16067>:	addps  %xmm12,%xmm13
   0x0000000000003ee7 <+16103>:	movaps 0xe50(%r8),%xmm6
   0x0000000000003eef <+16111>:	mulps  %xmm9,%xmm6
   0x0000000000003f03 <+16131>:	addps  %xmm12,%xmm6
   0x0000000000003f27 <+16167>:	movaps 0xe90(%r8),%xmm0
   0x0000000000003f2f <+16175>:	mulps  %xmm9,%xmm0
   0x0000000000003f43 <+16195>:	addps  %xmm12,%xmm0
   0x0000000000003f73 <+16243>:	movaps 0xed0(%r8),%xmm11
   0x0000000000003f7b <+16251>:	mulps  %xmm9,%xmm11
   0x0000000000003f9f <+16287>:	addps  %xmm12,%xmm11

1365
1366	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1367	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
1368	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003ea3 <+16035>:	movaps 0xe20(%r8),%xmm8
   0x0000000000003eaf <+16047>:	movaps 0x220(%rcx),%xmm10
   0x0000000000003ebf <+16063>:	mulps  %xmm10,%xmm8
   0x0000000000003ed3 <+16083>:	addps  %xmm13,%xmm8
   0x0000000000003ef7 <+16119>:	movaps 0xe60(%r8),%xmm8
   0x0000000000003eff <+16127>:	mulps  %xmm10,%xmm8
   0x0000000000003f13 <+16147>:	addps  %xmm6,%xmm8
   0x0000000000003f37 <+16183>:	movaps 0xea0(%r8),%xmm8
   0x0000000000003f3f <+16191>:	mulps  %xmm10,%xmm8
   0x0000000000003f5f <+16223>:	addps  %xmm0,%xmm8
   0x0000000000003f7f <+16255>:	movaps 0xee0(%r8),%xmm9
   0x0000000000003f87 <+16263>:	mulps  %xmm10,%xmm9
   0x0000000000003fab <+16299>:	addps  %xmm11,%xmm9

1369
1370	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1371	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
1372	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003eb7 <+16055>:	movaps 0x230(%rcx),%xmm15
   0x0000000000003ed7 <+16087>:	movaps 0xe30(%r8),%xmm13
   0x0000000000003edf <+16095>:	mulps  %xmm15,%xmm13
   0x0000000000003ef3 <+16115>:	addps  %xmm8,%xmm13
   0x0000000000003f17 <+16151>:	movaps 0xe70(%r8),%xmm6
   0x0000000000003f1f <+16159>:	mulps  %xmm15,%xmm6
   0x0000000000003f33 <+16179>:	addps  %xmm8,%xmm6
   0x0000000000003f53 <+16211>:	movaps 0xef0(%r8),%xmm14
   0x0000000000003f5b <+16219>:	mulps  %xmm15,%xmm14
   0x0000000000003f63 <+16227>:	movaps 0xeb0(%r8),%xmm0
   0x0000000000003f6b <+16235>:	mulps  %xmm15,%xmm0
   0x0000000000003f93 <+16275>:	addps  %xmm8,%xmm0
   0x0000000000003fc3 <+16323>:	addps  %xmm9,%xmm14

1373	      }
1374
1375	      /* A(4,4)*B(4,1) = C(4,1). */
1376	      for (i = 0; i < 4; i++)
1377	      {
1378	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1379	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
1380	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003f8b <+16267>:	movaps 0xf00(%r8),%xmm10
   0x0000000000003faf <+16303>:	movaps 0x300(%rcx),%xmm15
   0x0000000000003fbf <+16319>:	mulps  %xmm15,%xmm10
   0x0000000000003fd3 <+16339>:	addps  %xmm13,%xmm10
   0x0000000000003ff7 <+16375>:	movaps 0xf40(%r8),%xmm8
   0x0000000000003fff <+16383>:	mulps  %xmm15,%xmm8
   0x0000000000004013 <+16403>:	addps  %xmm6,%xmm8
   0x0000000000004037 <+16439>:	movaps 0xf80(%r8),%xmm6
   0x000000000000403f <+16447>:	mulps  %xmm15,%xmm6
   0x0000000000004053 <+16467>:	addps  %xmm0,%xmm6
   0x0000000000004074 <+16500>:	movaps 0xfc0(%r8),%xmm0
   0x000000000000407c <+16508>:	mulps  %xmm15,%xmm0
   0x00000000000040cc <+16588>:	addps  %xmm14,%xmm0

1381
1382	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1383	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
1384	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003f97 <+16279>:	movaps 0xf10(%r8),%xmm8
   0x0000000000003fb7 <+16311>:	movaps 0x310(%rcx),%xmm11
   0x0000000000003fc7 <+16327>:	mulps  %xmm11,%xmm8
   0x0000000000003fe3 <+16355>:	addps  %xmm10,%xmm8
   0x0000000000004017 <+16407>:	movaps 0xf50(%r8),%xmm6
   0x000000000000401f <+16415>:	mulps  %xmm11,%xmm6
   0x0000000000004023 <+16419>:	addps  %xmm8,%xmm6
   0x0000000000004056 <+16470>:	movaps 0xf90(%r8),%xmm0
   0x000000000000405e <+16478>:	mulps  %xmm11,%xmm0
   0x0000000000004062 <+16482>:	addps  %xmm6,%xmm0
   0x00000000000040d0 <+16592>:	movaps 0xfd0(%r8),%xmm14
   0x00000000000040d8 <+16600>:	mulps  %xmm11,%xmm14
   0x00000000000040dc <+16604>:	addps  %xmm0,%xmm14

1385
1386	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1387	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
1388	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003fa3 <+16291>:	movaps 0xf20(%r8),%xmm12
   0x0000000000003fd7 <+16343>:	movaps 0x320(%rcx),%xmm13
   0x0000000000003fdf <+16351>:	mulps  %xmm13,%xmm12
   0x0000000000003ff3 <+16371>:	addps  %xmm8,%xmm12
   0x0000000000004007 <+16391>:	movaps 0xf60(%r8),%xmm12
   0x000000000000400f <+16399>:	mulps  %xmm13,%xmm12
   0x0000000000004033 <+16435>:	addps  %xmm6,%xmm12
   0x0000000000004065 <+16485>:	movaps 0xfa0(%r8),%xmm6
   0x000000000000406d <+16493>:	mulps  %xmm13,%xmm6
   0x0000000000004071 <+16497>:	addps  %xmm0,%xmm6
   0x0000000000004080 <+16512>:	movaps 0xfe0(%r8),%xmm15
   0x0000000000004088 <+16520>:	mulps  %xmm13,%xmm15
   0x00000000000040ec <+16620>:	addps  %xmm14,%xmm15

1389
1390	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1391	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
1392	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003fcb <+16331>:	movaps 0xf30(%r8),%xmm9
   0x0000000000003fe7 <+16359>:	movaps 0x330(%rcx),%xmm10
   0x0000000000003fef <+16367>:	mulps  %xmm10,%xmm9
   0x0000000000004003 <+16387>:	addps  %xmm12,%xmm9
   0x0000000000004027 <+16423>:	movaps 0xf70(%r8),%xmm8
   0x000000000000402f <+16431>:	mulps  %xmm10,%xmm8
   0x0000000000004043 <+16451>:	addps  %xmm12,%xmm8
   0x0000000000004047 <+16455>:	movaps 0xfb0(%r8),%xmm12
   0x000000000000404f <+16463>:	mulps  %xmm10,%xmm12
   0x000000000000408c <+16524>:	addps  %xmm6,%xmm12
   0x00000000000040e0 <+16608>:	movaps 0xff0(%r8),%xmm0
   0x00000000000040e8 <+16616>:	mulps  %xmm10,%xmm0
   0x00000000000040f0 <+16624>:	addps  %xmm15,%xmm0

1393	      }
1394
1395	      /* Store C(4,1) block. */
1396	      for (i = 0; i < 4; i++)
1397	      {
1398	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000004090 <+16528>:	mulps  %xmm2,%xmm9
   0x00000000000040a4 <+16548>:	mulps  %xmm2,%xmm8
   0x00000000000040b8 <+16568>:	mulps  %xmm2,%xmm12
   0x00000000000040f4 <+16628>:	mulps  %xmm2,%xmm0

1399	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_41]), C_row[i]);
   0x0000000000004094 <+16532>:	addps  0x300(%r9),%xmm9
   0x00000000000040a8 <+16552>:	addps  0x310(%r9),%xmm8
   0x00000000000040bc <+16572>:	addps  0x320(%r9),%xmm12
   0x00000000000040f7 <+16631>:	addps  0x330(%r9),%xmm0

1400	        _mm_store_ps(&C[i*4+C_OFFSET_41], C_row[i]);
   0x000000000000409c <+16540>:	movaps %xmm9,0x300(%r9)
   0x00000000000040b0 <+16560>:	movaps %xmm8,0x310(%r9)
   0x00000000000040c4 <+16580>:	movaps %xmm12,0x320(%r9)
   0x00000000000040ff <+16639>:	movaps %xmm0,0x330(%r9)

1401	      }
1402	    }
1403
1404	    /* Reset C(4,2) matrix accumulators */
1405	    C_row[0] = _mm_setzero_ps();
1406	    C_row[1] = _mm_setzero_ps();
1407	    C_row[2] = _mm_setzero_ps();
1408	    C_row[3] = _mm_setzero_ps();
1409
1410	    if (norm[12]*norm[17] >= tolerance &&
   0x0000000000004107 <+16647>:	mulss  %xmm7,%xmm5
   0x000000000000410b <+16651>:	comiss %xmm1,%xmm5
   0x000000000000410e <+16654>:	jb     0x45ff <stream_kernel+17919>

1411	        norm[13]*norm[21] >= tolerance &&
   0x0000000000004114 <+16660>:	movss  0x4c(%rdx,%rsi,1),%xmm0
   0x000000000000411a <+16666>:	mulss  0x6c(%rdx,%rsi,1),%xmm0
   0x0000000000004120 <+16672>:	comiss %xmm1,%xmm0
   0x0000000000004123 <+16675>:	jb     0x45ff <stream_kernel+17919>

1412	        norm[14]*norm[25] >= tolerance &&
   0x0000000000004129 <+16681>:	movss  0x50(%rdx,%rsi,1),%xmm0
   0x000000000000412f <+16687>:	mulss  0x7c(%rdx,%rsi,1),%xmm0
   0x0000000000004135 <+16693>:	comiss %xmm1,%xmm0
   0x0000000000004138 <+16696>:	jb     0x45ff <stream_kernel+17919>

1413	        norm[15]*norm[29] >= tolerance)
   0x000000000000413e <+16702>:	movss  0x54(%rdx,%rsi,1),%xmm0
   0x0000000000004144 <+16708>:	mulss  0x8c(%rdx,%rsi,1),%xmm0
   0x000000000000414d <+16717>:	comiss %xmm1,%xmm0
   0x0000000000004150 <+16720>:	jb     0x45ff <stream_kernel+17919>

1414	    {
1415	      /* A(4,1)*B(1,2) = C(4,2). */
1416	      for (i = 0; i < 4; i++)
1417	      {
1418	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1419	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
1420	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004156 <+16726>:	movaps 0x40(%rcx),%xmm13
   0x000000000000416a <+16746>:	movaps 0xc00(%r8),%xmm6
   0x0000000000004182 <+16770>:	movaps 0xc40(%r8),%xmm14
   0x00000000000041a2 <+16802>:	mulps  %xmm13,%xmm6
   0x00000000000041ce <+16846>:	movaps 0xc80(%r8),%xmm12
   0x00000000000041d6 <+16854>:	mulps  %xmm13,%xmm14
   0x000000000000420a <+16906>:	mulps  %xmm13,%xmm12
   0x0000000000004216 <+16918>:	movaps 0xcc0(%r8),%xmm12
   0x000000000000423d <+16957>:	mulps  %xmm13,%xmm12

1421
1422	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1423	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
1424	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000415b <+16731>:	movaps 0x50(%rcx),%xmm8
   0x0000000000004172 <+16754>:	movaps 0xc10(%r8),%xmm9
   0x000000000000418a <+16778>:	movaps 0xc50(%r8),%xmm0
   0x00000000000041a6 <+16806>:	mulps  %xmm8,%xmm9
   0x00000000000041aa <+16810>:	addps  %xmm6,%xmm9
   0x00000000000041da <+16858>:	mulps  %xmm8,%xmm0
   0x00000000000041de <+16862>:	addps  %xmm14,%xmm0
   0x00000000000041e2 <+16866>:	movaps 0xc90(%r8),%xmm14
   0x000000000000420e <+16910>:	mulps  %xmm8,%xmm14
   0x0000000000004212 <+16914>:	addps  %xmm12,%xmm14
   0x0000000000004241 <+16961>:	movaps 0xcd0(%r8),%xmm13
   0x0000000000004249 <+16969>:	mulps  %xmm8,%xmm13
   0x0000000000004255 <+16981>:	addps  %xmm12,%xmm13

1425
1426	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1427	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
1428	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004160 <+16736>:	movaps 0x60(%rcx),%xmm11
   0x000000000000417a <+16762>:	movaps 0xc20(%r8),%xmm12
   0x0000000000004192 <+16786>:	movaps 0xc60(%r8),%xmm15
   0x00000000000041b6 <+16822>:	mulps  %xmm11,%xmm12
   0x00000000000041ba <+16826>:	addps  %xmm9,%xmm12
   0x00000000000041ea <+16874>:	mulps  %xmm11,%xmm15
   0x00000000000041ee <+16878>:	addps  %xmm0,%xmm15
   0x00000000000041f2 <+16882>:	movaps 0xca0(%r8),%xmm0
   0x000000000000421e <+16926>:	mulps  %xmm11,%xmm0
   0x0000000000004222 <+16930>:	addps  %xmm14,%xmm0
   0x000000000000424d <+16973>:	movaps 0xce0(%r8),%xmm8
   0x0000000000004259 <+16985>:	mulps  %xmm11,%xmm8
   0x0000000000004265 <+16997>:	addps  %xmm13,%xmm8

1429
1430	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1431	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
1432	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004165 <+16741>:	movaps 0x70(%rcx),%xmm10
   0x000000000000419a <+16794>:	movaps 0xc70(%r8),%xmm5
   0x00000000000041ae <+16814>:	movaps 0xcb0(%r8),%xmm6
   0x00000000000041be <+16830>:	movaps 0xc30(%r8),%xmm9
   0x00000000000041c6 <+16838>:	mulps  %xmm10,%xmm9
   0x00000000000041ca <+16842>:	addps  %xmm12,%xmm9
   0x00000000000041fa <+16890>:	mulps  %xmm10,%xmm5
   0x00000000000041fe <+16894>:	addps  %xmm15,%xmm5
   0x000000000000422e <+16942>:	mulps  %xmm10,%xmm6
   0x0000000000004232 <+16946>:	addps  %xmm0,%xmm6
   0x0000000000004235 <+16949>:	movaps 0xcf0(%r8),%xmm0
   0x0000000000004271 <+17009>:	mulps  %xmm10,%xmm0
   0x0000000000004281 <+17025>:	addps  %xmm8,%xmm0

1433	      }
1434
1435	      /* A(4,2)*B(2,2) = C(4,2). */
1436	      for (i = 0; i < 4; i++)
1437	      {
1438	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1439	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
1440	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004226 <+16934>:	movaps 0x140(%rcx),%xmm14
   0x000000000000425d <+16989>:	movaps 0xd00(%r8),%xmm11
   0x000000000000427d <+17021>:	mulps  %xmm14,%xmm11
   0x0000000000004291 <+17041>:	addps  %xmm9,%xmm11
   0x00000000000042b1 <+17073>:	movaps 0xd40(%r8),%xmm11
   0x00000000000042c5 <+17093>:	mulps  %xmm14,%xmm11
   0x00000000000042cd <+17101>:	movaps 0xd80(%r8),%xmm13
   0x00000000000042d9 <+17113>:	addps  %xmm5,%xmm11
   0x00000000000042fd <+17149>:	mulps  %xmm14,%xmm13
   0x0000000000004301 <+17153>:	movaps 0xdc0(%r8),%xmm5
   0x0000000000004319 <+17177>:	addps  %xmm6,%xmm13
   0x000000000000433d <+17213>:	mulps  %xmm14,%xmm5
   0x0000000000004359 <+17241>:	addps  %xmm0,%xmm5

1441
1442	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1443	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
1444	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004275 <+17013>:	movaps 0xd10(%r8),%xmm10
   0x0000000000004295 <+17045>:	movaps 0x150(%rcx),%xmm9
   0x00000000000042a5 <+17061>:	mulps  %xmm9,%xmm10
   0x00000000000042a9 <+17065>:	addps  %xmm11,%xmm10
   0x00000000000042dd <+17117>:	movaps 0xd50(%r8),%xmm5
   0x00000000000042e5 <+17125>:	mulps  %xmm9,%xmm5
   0x00000000000042e9 <+17129>:	addps  %xmm11,%xmm5
   0x000000000000431d <+17181>:	movaps 0xd90(%r8),%xmm6
   0x0000000000004325 <+17189>:	mulps  %xmm9,%xmm6
   0x0000000000004329 <+17193>:	addps  %xmm13,%xmm6
   0x0000000000004341 <+17217>:	movaps 0xdd0(%r8),%xmm14
   0x000000000000434d <+17229>:	mulps  %xmm9,%xmm14
   0x000000000000437c <+17276>:	addps  %xmm5,%xmm14

1445
1446	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1447	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
1448	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004202 <+16898>:	movaps 0x160(%rcx),%xmm15
   0x0000000000004269 <+17001>:	movaps 0xd20(%r8),%xmm13
   0x000000000000428d <+17037>:	mulps  %xmm15,%xmm13
   0x00000000000042b9 <+17081>:	addps  %xmm10,%xmm13
   0x00000000000042ed <+17133>:	movaps 0xd60(%r8),%xmm11
   0x00000000000042f5 <+17141>:	mulps  %xmm15,%xmm11
   0x00000000000042f9 <+17145>:	addps  %xmm5,%xmm11
   0x000000000000432d <+17197>:	movaps 0xda0(%r8),%xmm13
   0x0000000000004335 <+17205>:	mulps  %xmm15,%xmm13
   0x0000000000004339 <+17209>:	addps  %xmm6,%xmm13
   0x000000000000435c <+17244>:	movaps 0xde0(%r8),%xmm0
   0x0000000000004364 <+17252>:	mulps  %xmm15,%xmm0
   0x000000000000439a <+17306>:	addps  %xmm14,%xmm0

1449
1450	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1451	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
1452	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004285 <+17029>:	movaps 0xd30(%r8),%xmm8
   0x000000000000429d <+17053>:	movaps 0x170(%rcx),%xmm12
   0x00000000000042ad <+17069>:	mulps  %xmm12,%xmm8
   0x00000000000042bd <+17085>:	movaps 0xd70(%r8),%xmm10
   0x00000000000042c9 <+17097>:	addps  %xmm13,%xmm8
   0x00000000000042d5 <+17109>:	mulps  %xmm12,%xmm10
   0x0000000000004309 <+17161>:	addps  %xmm11,%xmm10
   0x000000000000430d <+17165>:	movaps 0xdb0(%r8),%xmm11
   0x0000000000004315 <+17173>:	mulps  %xmm12,%xmm11
   0x0000000000004349 <+17225>:	addps  %xmm13,%xmm11
   0x0000000000004351 <+17233>:	movaps 0xdf0(%r8),%xmm9
   0x0000000000004370 <+17264>:	mulps  %xmm12,%xmm9
   0x00000000000043aa <+17322>:	addps  %xmm0,%xmm9

1453	      }
1454
1455	      /* A(4,3)*B(3,2) = C(4,2). */
1456	      for (i = 0; i < 4; i++)
1457	      {
1458	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1459	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
1460	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004374 <+17268>:	movaps 0xe00(%r8),%xmm12
   0x00000000000043ae <+17326>:	movaps 0x240(%rcx),%xmm0
   0x00000000000043b5 <+17333>:	mulps  %xmm0,%xmm12
   0x00000000000043b9 <+17337>:	addps  %xmm8,%xmm12
   0x00000000000043bd <+17341>:	movaps 0xe40(%r8),%xmm8
   0x00000000000043c5 <+17349>:	mulps  %xmm0,%xmm8
   0x00000000000043dd <+17373>:	movaps 0xe80(%r8),%xmm15
   0x00000000000043e5 <+17381>:	mulps  %xmm0,%xmm15
   0x00000000000043f9 <+17401>:	addps  %xmm10,%xmm8
   0x000000000000440d <+17421>:	movaps 0xec0(%r8),%xmm8
   0x0000000000004415 <+17429>:	mulps  %xmm0,%xmm8
   0x0000000000004429 <+17449>:	addps  %xmm11,%xmm15
   0x0000000000004489 <+17545>:	addps  %xmm9,%xmm8

1461
1462	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1463	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
1464	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004368 <+17256>:	movaps 0xe10(%r8),%xmm15
   0x0000000000004380 <+17280>:	movaps 0x250(%rcx),%xmm5
   0x0000000000004396 <+17302>:	mulps  %xmm5,%xmm15
   0x00000000000043c9 <+17353>:	addps  %xmm12,%xmm15
   0x00000000000043ed <+17389>:	movaps 0xe50(%r8),%xmm14
   0x00000000000043f5 <+17397>:	mulps  %xmm5,%xmm14
   0x0000000000004409 <+17417>:	addps  %xmm8,%xmm14
   0x000000000000442d <+17453>:	movaps 0xe90(%r8),%xmm11
   0x0000000000004435 <+17461>:	mulps  %xmm5,%xmm11
   0x0000000000004449 <+17481>:	addps  %xmm15,%xmm11
   0x000000000000448d <+17549>:	movaps 0xed0(%r8),%xmm9
   0x0000000000004495 <+17557>:	mulps  %xmm5,%xmm9
   0x0000000000004499 <+17561>:	addps  %xmm8,%xmm9

1465
1466	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1467	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
1468	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004387 <+17287>:	movaps 0x260(%rcx),%xmm13
   0x000000000000439e <+17310>:	movaps 0xe20(%r8),%xmm14
   0x00000000000043a6 <+17318>:	mulps  %xmm13,%xmm14
   0x00000000000043d9 <+17369>:	addps  %xmm15,%xmm14
   0x00000000000043fd <+17405>:	movaps 0xe60(%r8),%xmm10
   0x0000000000004405 <+17413>:	mulps  %xmm13,%xmm10
   0x0000000000004419 <+17433>:	addps  %xmm14,%xmm10
   0x000000000000443d <+17469>:	movaps 0xea0(%r8),%xmm10
   0x0000000000004445 <+17477>:	mulps  %xmm13,%xmm10
   0x0000000000004459 <+17497>:	addps  %xmm11,%xmm10
   0x000000000000445d <+17501>:	movaps 0xee0(%r8),%xmm11
   0x0000000000004465 <+17509>:	mulps  %xmm13,%xmm11
   0x00000000000044b0 <+17584>:	addps  %xmm9,%xmm11

1469
1470	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1471	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
1472	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000438f <+17295>:	movaps 0x270(%rcx),%xmm6
   0x00000000000043cd <+17357>:	movaps 0xe30(%r8),%xmm12
   0x00000000000043d5 <+17365>:	mulps  %xmm6,%xmm12
   0x00000000000043e9 <+17385>:	addps  %xmm14,%xmm12
   0x000000000000441d <+17437>:	movaps 0xe70(%r8),%xmm14
   0x0000000000004425 <+17445>:	mulps  %xmm6,%xmm14
   0x0000000000004439 <+17465>:	addps  %xmm10,%xmm14
   0x000000000000444d <+17485>:	movaps 0xeb0(%r8),%xmm15
   0x0000000000004455 <+17493>:	mulps  %xmm6,%xmm15
   0x0000000000004469 <+17513>:	movaps 0xef0(%r8),%xmm13
   0x0000000000004471 <+17521>:	mulps  %xmm6,%xmm13
   0x0000000000004475 <+17525>:	addps  %xmm10,%xmm15
   0x00000000000044bf <+17599>:	addps  %xmm11,%xmm13

1473	      }
1474
1475	      /* A(4,4)*B(4,2) = C(4,2). */
1476	      for (i = 0; i < 4; i++)
1477	      {
1478	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1479	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
1480	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004479 <+17529>:	movaps 0xf00(%r8),%xmm10
   0x0000000000004481 <+17537>:	movaps 0xf40(%r8),%xmm6
   0x000000000000449d <+17565>:	movaps 0x340(%rcx),%xmm0
   0x00000000000044ac <+17580>:	mulps  %xmm0,%xmm10
   0x00000000000044bc <+17596>:	mulps  %xmm0,%xmm6
   0x00000000000044d7 <+17623>:	addps  %xmm12,%xmm10
   0x00000000000044e7 <+17639>:	addps  %xmm14,%xmm6
   0x0000000000004517 <+17687>:	movaps 0xf80(%r8),%xmm5
   0x0000000000004527 <+17703>:	movaps 0xfc0(%r8),%xmm6
   0x000000000000452f <+17711>:	mulps  %xmm0,%xmm5
   0x0000000000004536 <+17718>:	mulps  %xmm0,%xmm6
   0x0000000000004541 <+17729>:	addps  %xmm15,%xmm5
   0x0000000000004559 <+17753>:	addps  %xmm13,%xmm6

1481
1482	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1483	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
1484	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000044a4 <+17572>:	movaps 0x350(%rcx),%xmm8
   0x00000000000044db <+17627>:	movaps 0xf10(%r8),%xmm12
   0x00000000000044e3 <+17635>:	mulps  %xmm8,%xmm12
   0x00000000000044f3 <+17651>:	addps  %xmm10,%xmm12
   0x000000000000450b <+17675>:	movaps 0xf50(%r8),%xmm12
   0x000000000000451f <+17695>:	mulps  %xmm8,%xmm12
   0x0000000000004523 <+17699>:	addps  %xmm6,%xmm12
   0x0000000000004539 <+17721>:	movaps 0xfd0(%r8),%xmm0
   0x0000000000004545 <+17733>:	movaps 0xf90(%r8),%xmm15
   0x000000000000455d <+17757>:	mulps  %xmm8,%xmm15
   0x0000000000004579 <+17785>:	addps  %xmm5,%xmm15
   0x000000000000458d <+17805>:	mulps  %xmm8,%xmm0
   0x0000000000004599 <+17817>:	addps  %xmm6,%xmm0

1485
1486	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1487	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
1488	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000044c3 <+17603>:	movaps 0xf20(%r8),%xmm5
   0x00000000000044eb <+17643>:	movaps 0xf60(%r8),%xmm14
   0x00000000000044f7 <+17655>:	movaps 0x360(%rcx),%xmm10
   0x00000000000044ff <+17663>:	mulps  %xmm10,%xmm5
   0x0000000000004503 <+17667>:	addps  %xmm12,%xmm5
   0x0000000000004507 <+17671>:	mulps  %xmm10,%xmm14
   0x0000000000004532 <+17714>:	addps  %xmm12,%xmm14
   0x0000000000004561 <+17761>:	movaps 0xfe0(%r8),%xmm13
   0x000000000000457d <+17789>:	movaps 0xfa0(%r8),%xmm5
   0x0000000000004585 <+17797>:	mulps  %xmm10,%xmm5
   0x0000000000004589 <+17801>:	addps  %xmm15,%xmm5
   0x0000000000004595 <+17813>:	mulps  %xmm10,%xmm13
   0x00000000000045d8 <+17880>:	addps  %xmm0,%xmm13

1489
1490	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1491	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
1492	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000044b4 <+17588>:	movaps 0x370(%rcx),%xmm9
   0x00000000000044cb <+17611>:	movaps 0xf30(%r8),%xmm11
   0x00000000000044d3 <+17619>:	mulps  %xmm9,%xmm11
   0x0000000000004513 <+17683>:	addps  %xmm5,%xmm11
   0x000000000000454d <+17741>:	movaps 0xf70(%r8),%xmm12
   0x0000000000004555 <+17749>:	mulps  %xmm9,%xmm12
   0x0000000000004569 <+17769>:	addps  %xmm14,%xmm12
   0x000000000000456d <+17773>:	movaps 0xfb0(%r8),%xmm14
   0x0000000000004575 <+17781>:	mulps  %xmm9,%xmm14
   0x0000000000004591 <+17809>:	addps  %xmm5,%xmm14
   0x00000000000045dc <+17884>:	movaps 0xff0(%r8),%xmm0
   0x00000000000045e4 <+17892>:	mulps  %xmm9,%xmm0
   0x00000000000045e8 <+17896>:	addps  %xmm13,%xmm0

1493	      }
1494
1495	      /* Store C(4,2) block. */
1496	      for (i = 0; i < 4; i++)
1497	      {
1498	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x000000000000459c <+17820>:	mulps  %xmm2,%xmm11
   0x00000000000045b0 <+17840>:	mulps  %xmm2,%xmm12
   0x00000000000045c4 <+17860>:	mulps  %xmm2,%xmm14
   0x00000000000045ec <+17900>:	mulps  %xmm2,%xmm0

1499	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_42]), C_row[i]);
   0x00000000000045a0 <+17824>:	addps  0x340(%r9),%xmm11
   0x00000000000045b4 <+17844>:	addps  0x350(%r9),%xmm12
   0x00000000000045c8 <+17864>:	addps  0x360(%r9),%xmm14
   0x00000000000045ef <+17903>:	addps  0x370(%r9),%xmm0

1500	        _mm_store_ps(&C[i*4+C_OFFSET_42], C_row[i]);
   0x00000000000045a8 <+17832>:	movaps %xmm11,0x340(%r9)
   0x00000000000045bc <+17852>:	movaps %xmm12,0x350(%r9)
   0x00000000000045d0 <+17872>:	movaps %xmm14,0x360(%r9)
   0x00000000000045f7 <+17911>:	movaps %xmm0,0x370(%r9)

1501	      }
1502	    }
1503
1504	    /* Reset C(4,3) matrix accumulators */
1505	    C_row[0] = _mm_setzero_ps();
1506	    C_row[1] = _mm_setzero_ps();
1507	    C_row[2] = _mm_setzero_ps();
1508	    C_row[3] = _mm_setzero_ps();
1509
1510	    if (norm[12]*norm[18] >= tolerance &&
   0x00000000000045ff <+17919>:	mulss  %xmm7,%xmm4
   0x0000000000004603 <+17923>:	comiss %xmm1,%xmm4
   0x0000000000004606 <+17926>:	jb     0x4aff <stream_kernel+19199>

1511	        norm[13]*norm[22] >= tolerance &&
   0x000000000000460c <+17932>:	movss  0x4c(%rdx,%rsi,1),%xmm0
   0x0000000000004612 <+17938>:	mulss  0x70(%rdx,%rsi,1),%xmm0
   0x0000000000004618 <+17944>:	comiss %xmm1,%xmm0
   0x000000000000461b <+17947>:	jb     0x4aff <stream_kernel+19199>

1512	        norm[14]*norm[26] >= tolerance &&
   0x0000000000004621 <+17953>:	movss  0x50(%rdx,%rsi,1),%xmm0
   0x0000000000004627 <+17959>:	mulss  0x80(%rdx,%rsi,1),%xmm0
   0x0000000000004630 <+17968>:	comiss %xmm1,%xmm0
   0x0000000000004633 <+17971>:	jb     0x4aff <stream_kernel+19199>

1513	        norm[15]*norm[30] >= tolerance)
   0x0000000000004639 <+17977>:	movss  0x54(%rdx,%rsi,1),%xmm0
   0x000000000000463f <+17983>:	mulss  0x90(%rdx,%rsi,1),%xmm0
   0x0000000000004648 <+17992>:	comiss %xmm1,%xmm0
   0x000000000000464b <+17995>:	jb     0x4aff <stream_kernel+19199>

1514	    {
1515	      /* A(4,1)*B(1,3) = C(4,3). */
1516	      for (i = 0; i < 4; i++)
1517	      {
1518	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1519	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
1520	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004651 <+18001>:	movaps 0x80(%rcx),%xmm5
   0x0000000000004670 <+18032>:	movaps 0xc00(%r8),%xmm0
   0x0000000000004690 <+18064>:	movaps 0xc40(%r8),%xmm11
   0x00000000000046a8 <+18088>:	mulps  %xmm5,%xmm0
   0x00000000000046af <+18095>:	movaps 0xcc0(%r8),%xmm8
   0x00000000000046cb <+18123>:	movaps 0xc80(%r8),%xmm12
   0x00000000000046e3 <+18147>:	mulps  %xmm5,%xmm11
   0x0000000000004716 <+18198>:	mulps  %xmm5,%xmm12
   0x000000000000474a <+18250>:	mulps  %xmm5,%xmm8

1521
1522	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1523	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
1524	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004658 <+18008>:	movaps 0x90(%rcx),%xmm15
   0x0000000000004678 <+18040>:	movaps 0xc10(%r8),%xmm12
   0x0000000000004698 <+18072>:	movaps 0xc50(%r8),%xmm9
   0x00000000000046ab <+18091>:	mulps  %xmm15,%xmm12
   0x00000000000046b7 <+18103>:	addps  %xmm0,%xmm12
   0x00000000000046db <+18139>:	movaps 0xc90(%r8),%xmm13
   0x00000000000046e7 <+18151>:	mulps  %xmm15,%xmm9
   0x00000000000046eb <+18155>:	addps  %xmm11,%xmm9
   0x000000000000471a <+18202>:	mulps  %xmm15,%xmm13
   0x000000000000471e <+18206>:	addps  %xmm12,%xmm13
   0x000000000000474e <+18254>:	movaps 0xcd0(%r8),%xmm5
   0x0000000000004756 <+18262>:	mulps  %xmm15,%xmm5
   0x0000000000004762 <+18274>:	addps  %xmm8,%xmm5

1525
1526	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1527	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
1528	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004660 <+18016>:	movaps 0xa0(%rcx),%xmm10
   0x0000000000004680 <+18048>:	movaps 0xc20(%r8),%xmm13
   0x00000000000046a0 <+18080>:	movaps 0xc60(%r8),%xmm4
   0x00000000000046bb <+18107>:	mulps  %xmm10,%xmm13
   0x00000000000046c7 <+18119>:	addps  %xmm12,%xmm13
   0x00000000000046ef <+18159>:	movaps 0xca0(%r8),%xmm11
   0x00000000000046f7 <+18167>:	mulps  %xmm10,%xmm4
   0x00000000000046fb <+18171>:	addps  %xmm9,%xmm4
   0x000000000000472a <+18218>:	mulps  %xmm10,%xmm11
   0x000000000000472e <+18222>:	addps  %xmm13,%xmm11
   0x000000000000475a <+18266>:	movaps 0xce0(%r8),%xmm15
   0x0000000000004766 <+18278>:	mulps  %xmm10,%xmm15
   0x0000000000004772 <+18290>:	addps  %xmm5,%xmm15

1529
1530	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1531	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
1532	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004668 <+18024>:	movaps 0xb0(%rcx),%xmm14
   0x0000000000004688 <+18056>:	movaps 0xc30(%r8),%xmm6
   0x00000000000046bf <+18111>:	movaps 0xc70(%r8),%xmm0
   0x00000000000046d3 <+18131>:	mulps  %xmm14,%xmm6
   0x00000000000046d7 <+18135>:	addps  %xmm13,%xmm6
   0x0000000000004707 <+18183>:	mulps  %xmm14,%xmm0
   0x000000000000470b <+18187>:	addps  %xmm4,%xmm0
   0x000000000000470e <+18190>:	movaps 0xcb0(%r8),%xmm4
   0x000000000000473a <+18234>:	mulps  %xmm14,%xmm4
   0x000000000000473e <+18238>:	addps  %xmm11,%xmm4
   0x0000000000004742 <+18242>:	movaps 0xcf0(%r8),%xmm11
   0x000000000000477e <+18302>:	mulps  %xmm14,%xmm11
   0x000000000000478e <+18318>:	addps  %xmm15,%xmm11

1533	      }
1534
1535	      /* A(4,2)*B(2,3) = C(4,3). */
1536	      for (i = 0; i < 4; i++)
1537	      {
1538	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1539	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
1540	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004732 <+18226>:	movaps 0x180(%rcx),%xmm13
   0x000000000000476a <+18282>:	movaps 0xd00(%r8),%xmm10
   0x0000000000004776 <+18294>:	movaps 0xd40(%r8),%xmm5
   0x000000000000478a <+18314>:	mulps  %xmm13,%xmm10
   0x00000000000047aa <+18346>:	addps  %xmm6,%xmm10
   0x00000000000047ae <+18350>:	mulps  %xmm13,%xmm5
   0x00000000000047ca <+18378>:	addps  %xmm0,%xmm5
   0x00000000000047dd <+18397>:	movaps 0xd80(%r8),%xmm14
   0x00000000000047e5 <+18405>:	mulps  %xmm13,%xmm14
   0x000000000000481c <+18460>:	movaps 0xdc0(%r8),%xmm10
   0x0000000000004828 <+18472>:	addps  %xmm4,%xmm14
   0x000000000000483c <+18492>:	mulps  %xmm13,%xmm10
   0x0000000000004876 <+18550>:	addps  %xmm11,%xmm10

1541
1542	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1543	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
1544	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004782 <+18306>:	movaps 0xd10(%r8),%xmm14
   0x0000000000004796 <+18326>:	movaps 0x190(%rcx),%xmm8
   0x00000000000047a6 <+18342>:	mulps  %xmm8,%xmm14
   0x00000000000047b2 <+18354>:	movaps 0xd50(%r8),%xmm6
   0x00000000000047ba <+18362>:	addps  %xmm10,%xmm14
   0x00000000000047cd <+18381>:	mulps  %xmm8,%xmm6
   0x00000000000047f9 <+18425>:	addps  %xmm5,%xmm6
   0x000000000000482c <+18476>:	movaps 0xd90(%r8),%xmm4
   0x0000000000004834 <+18484>:	mulps  %xmm8,%xmm4
   0x0000000000004838 <+18488>:	addps  %xmm14,%xmm4
   0x0000000000004840 <+18496>:	movaps 0xdd0(%r8),%xmm13
   0x0000000000004853 <+18515>:	mulps  %xmm8,%xmm13
   0x0000000000004892 <+18578>:	addps  %xmm10,%xmm13

1545
1546	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1547	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
1548	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000479e <+18334>:	movaps 0x1a0(%rcx),%xmm15
   0x00000000000047be <+18366>:	movaps 0xd20(%r8),%xmm10
   0x00000000000047c6 <+18374>:	mulps  %xmm15,%xmm10
   0x00000000000047d1 <+18385>:	movaps 0xda0(%r8),%xmm0
   0x00000000000047d9 <+18393>:	addps  %xmm14,%xmm10
   0x00000000000047ed <+18413>:	movaps 0xd60(%r8),%xmm10
   0x00000000000047f5 <+18421>:	mulps  %xmm15,%xmm10
   0x0000000000004808 <+18440>:	addps  %xmm6,%xmm10
   0x0000000000004814 <+18452>:	mulps  %xmm15,%xmm0
   0x0000000000004848 <+18504>:	addps  %xmm4,%xmm0
   0x000000000000487a <+18554>:	movaps 0xde0(%r8),%xmm11
   0x0000000000004882 <+18562>:	mulps  %xmm15,%xmm11
   0x00000000000048a6 <+18598>:	addps  %xmm13,%xmm11

1549
1550	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1551	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
1552	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000046ff <+18175>:	movaps 0x1b0(%rcx),%xmm9
   0x0000000000004722 <+18210>:	movaps 0xd30(%r8),%xmm12
   0x0000000000004792 <+18322>:	mulps  %xmm9,%xmm12
   0x00000000000047e9 <+18409>:	addps  %xmm10,%xmm12
   0x00000000000047fc <+18428>:	movaps 0xd70(%r8),%xmm5
   0x0000000000004804 <+18436>:	mulps  %xmm9,%xmm5
   0x000000000000480c <+18444>:	movaps 0xdb0(%r8),%xmm6
   0x0000000000004818 <+18456>:	addps  %xmm10,%xmm5
   0x0000000000004824 <+18468>:	mulps  %xmm9,%xmm6
   0x000000000000484b <+18507>:	movaps 0xdf0(%r8),%xmm4
   0x000000000000485f <+18527>:	mulps  %xmm9,%xmm4
   0x000000000000486b <+18539>:	addps  %xmm0,%xmm6
   0x00000000000048b6 <+18614>:	addps  %xmm11,%xmm4

1553	      }
1554
1555	      /* A(4,3)*B(3,3) = C(4,3). */
1556	      for (i = 0; i < 4; i++)
1557	      {
1558	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1559	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
1560	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000486e <+18542>:	movaps 0x280(%rcx),%xmm14
   0x0000000000004886 <+18566>:	movaps 0xe00(%r8),%xmm15
   0x000000000000488e <+18574>:	mulps  %xmm14,%xmm15
   0x00000000000048c6 <+18630>:	addps  %xmm12,%xmm15
   0x00000000000048de <+18654>:	movaps 0xe40(%r8),%xmm15
   0x00000000000048ea <+18666>:	mulps  %xmm14,%xmm15
   0x0000000000004906 <+18694>:	addps  %xmm5,%xmm15
   0x000000000000491a <+18714>:	movaps 0xe80(%r8),%xmm15
   0x0000000000004922 <+18722>:	mulps  %xmm14,%xmm15
   0x0000000000004929 <+18729>:	movaps 0xec0(%r8),%xmm5
   0x0000000000004939 <+18745>:	mulps  %xmm14,%xmm5
   0x0000000000004945 <+18757>:	addps  %xmm6,%xmm15
   0x0000000000004955 <+18773>:	addps  %xmm4,%xmm5

1561
1562	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1563	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
1564	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004863 <+18531>:	movaps 0xe10(%r8),%xmm9
   0x00000000000048ba <+18618>:	movaps 0x290(%rcx),%xmm11
   0x00000000000048c2 <+18626>:	mulps  %xmm11,%xmm9
   0x00000000000048d6 <+18646>:	addps  %xmm15,%xmm9
   0x000000000000490a <+18698>:	movaps 0xe50(%r8),%xmm5
   0x0000000000004912 <+18706>:	mulps  %xmm11,%xmm5
   0x0000000000004916 <+18710>:	addps  %xmm15,%xmm5
   0x000000000000493d <+18749>:	movaps 0xed0(%r8),%xmm14
   0x0000000000004949 <+18761>:	movaps 0xe90(%r8),%xmm6
   0x0000000000004951 <+18769>:	mulps  %xmm11,%xmm6
   0x0000000000004958 <+18776>:	mulps  %xmm11,%xmm14
   0x0000000000004964 <+18788>:	addps  %xmm15,%xmm6
   0x0000000000004974 <+18804>:	addps  %xmm5,%xmm14

1565
1566	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1567	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
1568	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004857 <+18519>:	movaps 0xe20(%r8),%xmm8
   0x000000000000489e <+18590>:	movaps 0xe60(%r8),%xmm0
   0x00000000000048ca <+18634>:	movaps 0x2a0(%rcx),%xmm12
   0x00000000000048d2 <+18642>:	mulps  %xmm12,%xmm8
   0x00000000000048da <+18650>:	mulps  %xmm12,%xmm0
   0x00000000000048e6 <+18662>:	addps  %xmm9,%xmm8
   0x0000000000004926 <+18726>:	addps  %xmm5,%xmm0
   0x000000000000495c <+18780>:	movaps 0xee0(%r8),%xmm11
   0x0000000000004968 <+18792>:	movaps 0xea0(%r8),%xmm15
   0x0000000000004970 <+18800>:	mulps  %xmm12,%xmm15
   0x0000000000004978 <+18808>:	mulps  %xmm12,%xmm11
   0x00000000000049a0 <+18848>:	addps  %xmm6,%xmm15
   0x00000000000049be <+18878>:	addps  %xmm14,%xmm11

1569
1570	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1571	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
1572	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004896 <+18582>:	movaps 0xe30(%r8),%xmm10
   0x00000000000048aa <+18602>:	movaps 0x2b0(%rcx),%xmm13
   0x00000000000048b2 <+18610>:	mulps  %xmm13,%xmm10
   0x00000000000048ee <+18670>:	movaps 0xeb0(%r8),%xmm9
   0x00000000000048f6 <+18678>:	addps  %xmm8,%xmm10
   0x00000000000048fa <+18682>:	movaps 0xe70(%r8),%xmm8
   0x0000000000004902 <+18690>:	mulps  %xmm13,%xmm8
   0x0000000000004931 <+18737>:	mulps  %xmm13,%xmm9
   0x0000000000004935 <+18741>:	addps  %xmm0,%xmm8
   0x000000000000497c <+18812>:	movaps 0xef0(%r8),%xmm5
   0x000000000000499c <+18844>:	mulps  %xmm13,%xmm5
   0x00000000000049af <+18863>:	addps  %xmm15,%xmm9
   0x00000000000049cd <+18893>:	addps  %xmm11,%xmm5

1573	      }
1574
1575	      /* A(4,4)*B(4,3) = C(4,3). */
1576	      for (i = 0; i < 4; i++)
1577	      {
1578	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1579	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
1580	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004984 <+18820>:	movaps 0xf00(%r8),%xmm12
   0x000000000000498c <+18828>:	movaps 0xf40(%r8),%xmm4
   0x0000000000004994 <+18836>:	movaps 0xf80(%r8),%xmm0
   0x00000000000049a4 <+18852>:	movaps 0x380(%rcx),%xmm6
   0x00000000000049ab <+18859>:	mulps  %xmm6,%xmm12
   0x00000000000049bb <+18875>:	mulps  %xmm6,%xmm4
   0x00000000000049ca <+18890>:	mulps  %xmm6,%xmm0
   0x00000000000049e5 <+18917>:	addps  %xmm10,%xmm12
   0x0000000000004a15 <+18965>:	addps  %xmm8,%xmm4
   0x0000000000004a69 <+19049>:	addps  %xmm9,%xmm0
   0x0000000000004a8d <+19085>:	movaps 0xfc0(%r8),%xmm0
   0x0000000000004a95 <+19093>:	mulps  %xmm6,%xmm0
   0x0000000000004ad4 <+19156>:	addps  %xmm5,%xmm0

1581
1582	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1583	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
1584	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000049b3 <+18867>:	movaps 0x390(%rcx),%xmm15
   0x00000000000049d1 <+18897>:	movaps 0xf10(%r8),%xmm13
   0x00000000000049d9 <+18905>:	movaps 0xf50(%r8),%xmm11
   0x00000000000049e1 <+18913>:	mulps  %xmm15,%xmm13
   0x00000000000049f5 <+18933>:	addps  %xmm12,%xmm13
   0x00000000000049f9 <+18937>:	mulps  %xmm15,%xmm11
   0x0000000000004a35 <+18997>:	addps  %xmm4,%xmm11
   0x0000000000004a45 <+19013>:	movaps 0xfd0(%r8),%xmm4
   0x0000000000004a55 <+19029>:	mulps  %xmm15,%xmm4
   0x0000000000004a6d <+19053>:	movaps 0xf90(%r8),%xmm9
   0x0000000000004a75 <+19061>:	mulps  %xmm15,%xmm9
   0x0000000000004a89 <+19081>:	addps  %xmm0,%xmm9
   0x0000000000004ae3 <+19171>:	addps  %xmm0,%xmm4

1585
1586	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1587	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
1588	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000049c2 <+18882>:	movaps 0x3a0(%rcx),%xmm14
   0x00000000000049e9 <+18921>:	movaps 0xf20(%r8),%xmm10
   0x00000000000049f1 <+18929>:	mulps  %xmm14,%xmm10
   0x0000000000004a05 <+18949>:	addps  %xmm13,%xmm10
   0x0000000000004a19 <+18969>:	movaps 0xfa0(%r8),%xmm8
   0x0000000000004a21 <+18977>:	mulps  %xmm14,%xmm8
   0x0000000000004a29 <+18985>:	movaps 0xf60(%r8),%xmm10
   0x0000000000004a31 <+18993>:	mulps  %xmm14,%xmm10
   0x0000000000004a59 <+19033>:	addps  %xmm11,%xmm10
   0x0000000000004aa0 <+19104>:	addps  %xmm9,%xmm8
   0x0000000000004ad7 <+19159>:	movaps 0xfe0(%r8),%xmm5
   0x0000000000004adf <+19167>:	mulps  %xmm14,%xmm5
   0x0000000000004ae6 <+19174>:	addps  %xmm4,%xmm5

1589
1590	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1591	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
1592	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000049fd <+18941>:	movaps 0x3b0(%rcx),%xmm12
   0x0000000000004a09 <+18953>:	movaps 0xf30(%r8),%xmm13
   0x0000000000004a11 <+18961>:	mulps  %xmm12,%xmm13
   0x0000000000004a25 <+18981>:	addps  %xmm10,%xmm13
   0x0000000000004a5d <+19037>:	movaps 0xf70(%r8),%xmm11
   0x0000000000004a65 <+19045>:	mulps  %xmm12,%xmm11
   0x0000000000004a79 <+19065>:	addps  %xmm10,%xmm11
   0x0000000000004a7d <+19069>:	movaps 0xfb0(%r8),%xmm10
   0x0000000000004a85 <+19077>:	mulps  %xmm12,%xmm10
   0x0000000000004a98 <+19096>:	movaps 0xff0(%r8),%xmm6
   0x0000000000004aa4 <+19108>:	mulps  %xmm12,%xmm6
   0x0000000000004aa8 <+19112>:	addps  %xmm8,%xmm10
   0x0000000000004ae9 <+19177>:	addps  %xmm5,%xmm6

1593	      }
1594
1595	      /* Store C(4,3) block. */
1596	      for (i = 0; i < 4; i++)
1597	      {
1598	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000004a39 <+19001>:	mulps  %xmm2,%xmm13
   0x0000000000004aac <+19116>:	mulps  %xmm2,%xmm11
   0x0000000000004ac0 <+19136>:	mulps  %xmm2,%xmm10
   0x0000000000004aec <+19180>:	mulps  %xmm2,%xmm6

1599	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_43]), C_row[i]);
   0x0000000000004a3d <+19005>:	addps  0x380(%r9),%xmm13
   0x0000000000004ab0 <+19120>:	addps  0x390(%r9),%xmm11
   0x0000000000004ac4 <+19140>:	addps  0x3a0(%r9),%xmm10
   0x0000000000004aef <+19183>:	addps  0x3b0(%r9),%xmm6

1600	        _mm_store_ps(&C[i*4+C_OFFSET_43], C_row[i]);
   0x0000000000004a4d <+19021>:	movaps %xmm13,0x380(%r9)
   0x0000000000004ab8 <+19128>:	movaps %xmm11,0x390(%r9)
   0x0000000000004acc <+19148>:	movaps %xmm10,0x3a0(%r9)
   0x0000000000004af7 <+19191>:	movaps %xmm6,0x3b0(%r9)

1601	      }
1602	    }
1603
1604	    /* Reset C(4,4) matrix accumulators */
1605	    C_row[0] = _mm_setzero_ps();
1606	    C_row[1] = _mm_setzero_ps();
1607	    C_row[2] = _mm_setzero_ps();
1608	    C_row[3] = _mm_setzero_ps();
1609
1610	    if (norm[12]*norm[19] >= tolerance &&
   0x0000000000004aff <+19199>:	mulss  %xmm7,%xmm3
   0x0000000000004b03 <+19203>:	comiss %xmm1,%xmm3
   0x0000000000004b06 <+19206>:	jb     0x4fff <stream_kernel+20479>

1611	        norm[13]*norm[23] >= tolerance &&
   0x0000000000004b0c <+19212>:	movss  0x4c(%rdx,%rsi,1),%xmm0
   0x0000000000004b12 <+19218>:	mulss  0x74(%rdx,%rsi,1),%xmm0
   0x0000000000004b18 <+19224>:	comiss %xmm1,%xmm0
   0x0000000000004b1b <+19227>:	jb     0x4fff <stream_kernel+20479>

1612	        norm[14]*norm[27] >= tolerance &&
   0x0000000000004b21 <+19233>:	movss  0x50(%rdx,%rsi,1),%xmm0
   0x0000000000004b27 <+19239>:	mulss  0x84(%rdx,%rsi,1),%xmm0
   0x0000000000004b30 <+19248>:	comiss %xmm1,%xmm0
   0x0000000000004b33 <+19251>:	jb     0x4fff <stream_kernel+20479>

1613	        norm[15]*norm[31] >= tolerance)
   0x0000000000004b39 <+19257>:	movss  0x54(%rdx,%rsi,1),%xmm0
   0x0000000000004b3f <+19263>:	mulss  0x94(%rdx,%rsi,1),%xmm0
   0x0000000000004b48 <+19272>:	comiss %xmm1,%xmm0
   0x0000000000004b4b <+19275>:	jb     0x4fff <stream_kernel+20479>

1614	    {
1615	      /* A(4,1)*B(1,4) = C(4,4). */
1616	      for (i = 0; i < 4; i++)
1617	      {
1618	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1619	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
1620	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004b51 <+19281>:	movaps 0xc0(%rcx),%xmm0
   0x0000000000004b70 <+19312>:	movaps 0xc00(%r8),%xmm6
   0x0000000000004b90 <+19344>:	movaps 0xc40(%r8),%xmm12
   0x0000000000004bb0 <+19376>:	movaps 0xc80(%r8),%xmm4
   0x0000000000004bc0 <+19392>:	mulps  %xmm0,%xmm6
   0x0000000000004bdb <+19419>:	movaps 0xcc0(%r8),%xmm9
   0x0000000000004bf3 <+19443>:	mulps  %xmm0,%xmm12
   0x0000000000004c26 <+19494>:	mulps  %xmm0,%xmm4
   0x0000000000004c56 <+19542>:	mulps  %xmm0,%xmm9

1621
1622	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1623	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
1624	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004b58 <+19288>:	movaps 0xd0(%rcx),%xmm15
   0x0000000000004b78 <+19320>:	movaps 0xc10(%r8),%xmm9
   0x0000000000004b98 <+19352>:	movaps 0xc50(%r8),%xmm11
   0x0000000000004bb8 <+19384>:	movaps 0xc90(%r8),%xmm3
   0x0000000000004bc3 <+19395>:	mulps  %xmm15,%xmm9
   0x0000000000004bc7 <+19399>:	addps  %xmm6,%xmm9
   0x0000000000004beb <+19435>:	movaps 0xcd0(%r8),%xmm10
   0x0000000000004bf7 <+19447>:	mulps  %xmm15,%xmm11
   0x0000000000004bfb <+19451>:	addps  %xmm12,%xmm11
   0x0000000000004c29 <+19497>:	mulps  %xmm15,%xmm3
   0x0000000000004c2d <+19501>:	addps  %xmm4,%xmm3
   0x0000000000004c62 <+19554>:	mulps  %xmm15,%xmm10
   0x0000000000004c6e <+19566>:	addps  %xmm9,%xmm10

1625
1626	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1627	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
1628	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004b60 <+19296>:	movaps 0xe0(%rcx),%xmm13
   0x0000000000004b80 <+19328>:	movaps 0xc20(%r8),%xmm10
   0x0000000000004ba0 <+19360>:	movaps 0xc60(%r8),%xmm5
   0x0000000000004bcb <+19403>:	movaps 0xca0(%r8),%xmm6
   0x0000000000004bd3 <+19411>:	mulps  %xmm13,%xmm10
   0x0000000000004bd7 <+19415>:	addps  %xmm9,%xmm10
   0x0000000000004c07 <+19463>:	mulps  %xmm13,%xmm5
   0x0000000000004c0b <+19467>:	addps  %xmm11,%xmm5
   0x0000000000004c38 <+19512>:	mulps  %xmm13,%xmm6
   0x0000000000004c3c <+19516>:	addps  %xmm3,%xmm6
   0x0000000000004c66 <+19558>:	movaps 0xce0(%r8),%xmm15
   0x0000000000004c72 <+19570>:	mulps  %xmm13,%xmm15
   0x0000000000004c7e <+19582>:	addps  %xmm10,%xmm15

1629
1630	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1631	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
1632	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004b68 <+19304>:	movaps 0xf0(%rcx),%xmm14
   0x0000000000004b88 <+19336>:	movaps 0xc30(%r8),%xmm7
   0x0000000000004ba8 <+19368>:	movaps 0xc70(%r8),%xmm8
   0x0000000000004be3 <+19427>:	mulps  %xmm14,%xmm7
   0x0000000000004be7 <+19431>:	addps  %xmm10,%xmm7
   0x0000000000004c17 <+19479>:	mulps  %xmm14,%xmm8
   0x0000000000004c1b <+19483>:	addps  %xmm5,%xmm8
   0x0000000000004c30 <+19504>:	movaps 0xcf0(%r8),%xmm4
   0x0000000000004c3f <+19519>:	movaps 0xcb0(%r8),%xmm3
   0x0000000000004c47 <+19527>:	mulps  %xmm14,%xmm3
   0x0000000000004c4b <+19531>:	addps  %xmm6,%xmm3
   0x0000000000004c82 <+19586>:	mulps  %xmm14,%xmm4
   0x0000000000004c8e <+19598>:	addps  %xmm15,%xmm4

1633	      }
1634
1635	      /* A(4,2)*B(2,4) = C(4,4). */
1636	      for (i = 0; i < 4; i++)
1637	      {
1638	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1639	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
1640	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004c1f <+19487>:	movaps 0x1c0(%rcx),%xmm5
   0x0000000000004c4e <+19534>:	movaps 0xd40(%r8),%xmm6
   0x0000000000004c76 <+19574>:	movaps 0xd00(%r8),%xmm13
   0x0000000000004c92 <+19602>:	mulps  %xmm5,%xmm13
   0x0000000000004caa <+19626>:	addps  %xmm7,%xmm13
   0x0000000000004cc6 <+19654>:	movaps 0xd80(%r8),%xmm15
   0x0000000000004cd6 <+19670>:	mulps  %xmm5,%xmm6
   0x0000000000004cf1 <+19697>:	addps  %xmm8,%xmm6
   0x0000000000004d11 <+19729>:	mulps  %xmm5,%xmm15
   0x0000000000004d19 <+19737>:	movaps 0xdc0(%r8),%xmm8
   0x0000000000004d31 <+19761>:	mulps  %xmm5,%xmm8
   0x0000000000004d3d <+19773>:	addps  %xmm3,%xmm15
   0x0000000000004d4d <+19789>:	addps  %xmm4,%xmm8

1641
1642	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1643	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
1644	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004bff <+19455>:	movaps 0x1d0(%rcx),%xmm12
   0x0000000000004c86 <+19590>:	movaps 0xd10(%r8),%xmm14
   0x0000000000004ca6 <+19622>:	mulps  %xmm12,%xmm14
   0x0000000000004cba <+19642>:	addps  %xmm13,%xmm14
   0x0000000000004cf9 <+19705>:	movaps 0xd50(%r8),%xmm8
   0x0000000000004d01 <+19713>:	mulps  %xmm12,%xmm8
   0x0000000000004d05 <+19717>:	addps  %xmm6,%xmm8
   0x0000000000004d35 <+19765>:	movaps 0xdd0(%r8),%xmm5
   0x0000000000004d41 <+19777>:	movaps 0xd90(%r8),%xmm3
   0x0000000000004d49 <+19785>:	mulps  %xmm12,%xmm3
   0x0000000000004d51 <+19793>:	mulps  %xmm12,%xmm5
   0x0000000000004d5d <+19805>:	addps  %xmm15,%xmm3
   0x0000000000004d7c <+19836>:	addps  %xmm8,%xmm5

1645
1646	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1647	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
1648	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004c96 <+19606>:	movaps 0x1e0(%rcx),%xmm9
   0x0000000000004cb2 <+19634>:	movaps 0xd20(%r8),%xmm7
   0x0000000000004cbe <+19646>:	movaps 0xd60(%r8),%xmm13
   0x0000000000004cce <+19662>:	mulps  %xmm9,%xmm7
   0x0000000000004cd2 <+19666>:	addps  %xmm14,%xmm7
   0x0000000000004ce5 <+19685>:	mulps  %xmm9,%xmm13
   0x0000000000004ce9 <+19689>:	movaps 0xda0(%r8),%xmm7
   0x0000000000004d15 <+19733>:	addps  %xmm8,%xmm13
   0x0000000000004d21 <+19745>:	mulps  %xmm9,%xmm7
   0x0000000000004d55 <+19797>:	movaps 0xde0(%r8),%xmm12
   0x0000000000004d61 <+19809>:	mulps  %xmm9,%xmm12
   0x0000000000004d6d <+19821>:	addps  %xmm3,%xmm7
   0x0000000000004db7 <+19895>:	addps  %xmm5,%xmm12

1649
1650	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1651	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
1652	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004c0f <+19471>:	movaps 0xd30(%r8),%xmm11
   0x0000000000004c5a <+19546>:	movaps 0xd70(%r8),%xmm0
   0x0000000000004c9e <+19614>:	movaps 0x1f0(%rcx),%xmm10
   0x0000000000004cae <+19630>:	mulps  %xmm10,%xmm11
   0x0000000000004ce1 <+19681>:	addps  %xmm7,%xmm11
   0x0000000000004cf5 <+19701>:	mulps  %xmm10,%xmm0
   0x0000000000004d09 <+19721>:	movaps 0xdf0(%r8),%xmm6
   0x0000000000004d25 <+19749>:	addps  %xmm13,%xmm0
   0x0000000000004d70 <+19824>:	mulps  %xmm10,%xmm6
   0x0000000000004d74 <+19828>:	movaps 0xdb0(%r8),%xmm3
   0x0000000000004d80 <+19840>:	mulps  %xmm10,%xmm3
   0x0000000000004d98 <+19864>:	addps  %xmm7,%xmm3
   0x0000000000004dc7 <+19911>:	addps  %xmm12,%xmm6

1653	      }
1654
1655	      /* A(4,3)*B(3,4) = C(4,4). */
1656	      for (i = 0; i < 4; i++)
1657	      {
1658	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1659	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
1660	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004d65 <+19813>:	movaps 0xe00(%r8),%xmm9
   0x0000000000004d84 <+19844>:	movaps 0x2c0(%rcx),%xmm8
   0x0000000000004d94 <+19860>:	mulps  %xmm8,%xmm9
   0x0000000000004da3 <+19875>:	movaps 0xe40(%r8),%xmm10
   0x0000000000004dc3 <+19907>:	mulps  %xmm8,%xmm10
   0x0000000000004dd7 <+19927>:	addps  %xmm11,%xmm9
   0x0000000000004de7 <+19943>:	addps  %xmm0,%xmm10
   0x0000000000004e2b <+20011>:	movaps 0xe80(%r8),%xmm10
   0x0000000000004e33 <+20019>:	mulps  %xmm8,%xmm10
   0x0000000000004e57 <+20055>:	addps  %xmm3,%xmm10
   0x0000000000004e5f <+20063>:	movaps 0xec0(%r8),%xmm3
   0x0000000000004e6b <+20075>:	mulps  %xmm8,%xmm3
   0x0000000000004ea3 <+20131>:	addps  %xmm6,%xmm3

1661
1662	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1663	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
1664	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004d8c <+19852>:	movaps 0x2d0(%rcx),%xmm15
   0x0000000000004dcb <+19915>:	movaps 0xe90(%r8),%xmm12
   0x0000000000004ddb <+19931>:	mulps  %xmm15,%xmm12
   0x0000000000004ddf <+19935>:	movaps 0xe10(%r8),%xmm11
   0x0000000000004deb <+19947>:	movaps 0xe50(%r8),%xmm0
   0x0000000000004df3 <+19955>:	mulps  %xmm15,%xmm11
   0x0000000000004df7 <+19959>:	addps  %xmm9,%xmm11
   0x0000000000004e0b <+19979>:	mulps  %xmm15,%xmm0
   0x0000000000004e27 <+20007>:	addps  %xmm10,%xmm0
   0x0000000000004e67 <+20071>:	addps  %xmm10,%xmm12
   0x0000000000004e6f <+20079>:	movaps 0xed0(%r8),%xmm8
   0x0000000000004e83 <+20099>:	mulps  %xmm15,%xmm8
   0x0000000000004eba <+20154>:	addps  %xmm3,%xmm8

1665
1666	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1667	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
1668	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004d29 <+19753>:	movaps 0x2e0(%rcx),%xmm13
   0x0000000000004dfb <+19963>:	movaps 0xe20(%r8),%xmm9
   0x0000000000004e03 <+19971>:	mulps  %xmm13,%xmm9
   0x0000000000004e07 <+19975>:	addps  %xmm11,%xmm9
   0x0000000000004e0f <+19983>:	movaps 0xea0(%r8),%xmm11
   0x0000000000004e1b <+19995>:	movaps 0xe60(%r8),%xmm9
   0x0000000000004e23 <+20003>:	mulps  %xmm13,%xmm9
   0x0000000000004e37 <+20023>:	addps  %xmm0,%xmm9
   0x0000000000004e3b <+20027>:	mulps  %xmm13,%xmm11
   0x0000000000004e77 <+20087>:	addps  %xmm12,%xmm11
   0x0000000000004ea6 <+20134>:	movaps 0xee0(%r8),%xmm6
   0x0000000000004eae <+20142>:	mulps  %xmm13,%xmm6
   0x0000000000004ed2 <+20178>:	addps  %xmm8,%xmm6

1669
1670	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1671	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
1672	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004cd9 <+19673>:	movaps 0x2f0(%rcx),%xmm14
   0x0000000000004d9b <+19867>:	movaps 0xe30(%r8),%xmm7
   0x0000000000004dab <+19883>:	movaps 0xe70(%r8),%xmm4
   0x0000000000004db3 <+19891>:	mulps  %xmm14,%xmm7
   0x0000000000004dbb <+19899>:	movaps 0xeb0(%r8),%xmm5
   0x0000000000004dd3 <+19923>:	mulps  %xmm14,%xmm4
   0x0000000000004e17 <+19991>:	addps  %xmm9,%xmm7
   0x0000000000004e47 <+20039>:	addps  %xmm9,%xmm4
   0x0000000000004e4b <+20043>:	mulps  %xmm14,%xmm5
   0x0000000000004e7b <+20091>:	movaps 0xef0(%r8),%xmm12
   0x0000000000004e87 <+20103>:	addps  %xmm11,%xmm5
   0x0000000000004e8b <+20107>:	mulps  %xmm14,%xmm12
   0x0000000000004ee2 <+20194>:	addps  %xmm6,%xmm12

1673	      }
1674
1675	      /* A(4,4)*B(4,4) = C(4,4). */
1676	      for (i = 0; i < 4; i++)
1677	      {
1678	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1679	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
1680	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004e97 <+20119>:	movaps 0xf40(%r8),%xmm15
   0x0000000000004eb2 <+20146>:	movaps 0xf00(%r8),%xmm13
   0x0000000000004ebe <+20158>:	movaps 0x3c0(%rcx),%xmm11
   0x0000000000004ece <+20174>:	mulps  %xmm11,%xmm13
   0x0000000000004ed6 <+20182>:	mulps  %xmm11,%xmm15
   0x0000000000004eee <+20206>:	movaps 0xf80(%r8),%xmm3
   0x0000000000004ef6 <+20214>:	mulps  %xmm11,%xmm3
   0x0000000000004efa <+20218>:	addps  %xmm7,%xmm13
   0x0000000000004f0a <+20234>:	addps  %xmm4,%xmm15
   0x0000000000004f4c <+20300>:	movaps 0xfc0(%r8),%xmm14
   0x0000000000004f5b <+20315>:	mulps  %xmm11,%xmm14
   0x0000000000004f6e <+20334>:	addps  %xmm5,%xmm3
   0x0000000000004f85 <+20357>:	addps  %xmm12,%xmm14

1681
1682	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1683	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
1684	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004ec6 <+20166>:	movaps 0x3d0(%rcx),%xmm10
   0x0000000000004efe <+20222>:	movaps 0xf10(%r8),%xmm7
   0x0000000000004f06 <+20230>:	mulps  %xmm10,%xmm7
   0x0000000000004f0e <+20238>:	movaps 0xf50(%r8),%xmm4
   0x0000000000004f16 <+20246>:	mulps  %xmm10,%xmm4
   0x0000000000004f1a <+20250>:	addps  %xmm13,%xmm7
   0x0000000000004f39 <+20281>:	addps  %xmm15,%xmm4
   0x0000000000004f5f <+20319>:	movaps 0xfd0(%r8),%xmm11
   0x0000000000004f6a <+20330>:	mulps  %xmm10,%xmm11
   0x0000000000004f71 <+20337>:	movaps 0xf90(%r8),%xmm5
   0x0000000000004f79 <+20345>:	mulps  %xmm10,%xmm5
   0x0000000000004f95 <+20373>:	addps  %xmm3,%xmm5
   0x0000000000004fdf <+20447>:	addps  %xmm14,%xmm11

1685
1686	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1687	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
1688	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004e3f <+20031>:	movaps 0xf60(%r8),%xmm0
   0x0000000000004e4f <+20047>:	movaps 0x3e0(%rcx),%xmm9
   0x0000000000004e5b <+20059>:	mulps  %xmm9,%xmm0
   0x0000000000004e8f <+20111>:	movaps 0xf20(%r8),%xmm14
   0x0000000000004e9f <+20127>:	mulps  %xmm9,%xmm14
   0x0000000000004f1e <+20254>:	movaps 0xfa0(%r8),%xmm13
   0x0000000000004f26 <+20262>:	mulps  %xmm9,%xmm13
   0x0000000000004f2a <+20266>:	addps  %xmm7,%xmm14
   0x0000000000004f58 <+20312>:	addps  %xmm4,%xmm0
   0x0000000000004f7d <+20349>:	movaps 0xfe0(%r8),%xmm10
   0x0000000000004f89 <+20361>:	mulps  %xmm9,%xmm10
   0x0000000000004f9c <+20380>:	addps  %xmm5,%xmm13
   0x0000000000004fe3 <+20451>:	addps  %xmm11,%xmm10

1689
1690	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1691	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
1692	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004eda <+20186>:	movaps 0xf30(%r8),%xmm8
   0x0000000000004ee6 <+20198>:	movaps 0xf70(%r8),%xmm6
   0x0000000000004f2e <+20270>:	movaps 0x3f0(%rcx),%xmm7
   0x0000000000004f35 <+20277>:	mulps  %xmm7,%xmm8
   0x0000000000004f3d <+20285>:	mulps  %xmm7,%xmm6
   0x0000000000004f40 <+20288>:	movaps 0xfb0(%r8),%xmm15
   0x0000000000004f48 <+20296>:	addps  %xmm14,%xmm8
   0x0000000000004f54 <+20308>:	mulps  %xmm7,%xmm15
   0x0000000000004f67 <+20327>:	addps  %xmm0,%xmm6
   0x0000000000004f8d <+20365>:	movaps 0xff0(%r8),%xmm9
   0x0000000000004f98 <+20376>:	mulps  %xmm7,%xmm9
   0x0000000000004fc7 <+20423>:	addps  %xmm13,%xmm15
   0x0000000000004fe7 <+20455>:	addps  %xmm10,%xmm9

1693	      }
1694
1695	      /* Store C(4,4) block. */
1696	      for (i = 0; i < 4; i++)
1697	      {
1698	        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000004fa0 <+20384>:	mulps  %xmm2,%xmm8
   0x0000000000004fb4 <+20404>:	mulps  %xmm2,%xmm6
   0x0000000000004fcb <+20427>:	mulps  %xmm2,%xmm15
   0x0000000000004feb <+20459>:	mulps  %xmm2,%xmm9

1699	        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_44]), C_row[i]);
   0x0000000000004fa4 <+20388>:	addps  0x3c0(%r9),%xmm8
   0x0000000000004fb7 <+20407>:	addps  0x3d0(%r9),%xmm6
   0x0000000000004fcf <+20431>:	addps  0x3e0(%r9),%xmm15
   0x0000000000004fef <+20463>:	addps  0x3f0(%r9),%xmm9

1700	        _mm_store_ps(&C[i*4+C_OFFSET_44], C_row[i]);
   0x0000000000004fac <+20396>:	movaps %xmm8,0x3c0(%r9)
   0x0000000000004fbf <+20415>:	movaps %xmm6,0x3d0(%r9)
   0x0000000000004fd7 <+20439>:	movaps %xmm15,0x3e0(%r9)
   0x0000000000004ff7 <+20471>:	movaps %xmm9,0x3f0(%r9)

1701	      }
1702	    }
1703	  }
1704	}
   0x000000000000500d <+20493>:	retq
   0x000000000000500e <+20494>:	xchg   %ax,%ax

End of assembler dump.
