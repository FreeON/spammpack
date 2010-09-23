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
   0x0000000000000000 <+0>:	movaps %xmm0,%xmm6
   0x0000000000000003 <+3>:	shufps $0x0,%xmm6,%xmm6

95
96	  for (stream_index = 0; stream_index < max_stream_index; stream_index++)
   0x0000000000000007 <+7>:	xor    %edx,%edx
   0x000000000000000c <+12>:	test   %edi,%edi
   0x000000000000000e <+14>:	jbe    0x5593 <stream_kernel+21907>
   0x0000000000000014 <+20>:	xor    %eax,%eax
   0x0000000000000016 <+22>:	xorps  %xmm2,%xmm2
   0x0000000000005544 <+21828>:	inc    %eax
   0x0000000000005546 <+21830>:	mov    %eax,%edx
   0x0000000000005550 <+21840>:	mov    %edx,%eax
   0x0000000000005552 <+21842>:	cmp    %edi,%eax
   0x000000000000558d <+21901>:	jb     0x19 <stream_kernel+25>

97	  {
98	    /* Load pointers to matrix data blocks. */
99	    A = multiply_stream[stream_index].A_block;
   0x0000000000000019 <+25>:	imul   $0x98,%rdx,%rdx
   0x0000000000000020 <+32>:	mov    (%rsi,%rdx,1),%r8

100	    B = multiply_stream[stream_index].B_block;
   0x0000000000000024 <+36>:	mov    0x8(%rsi,%rdx,1),%rcx

101	    C = multiply_stream[stream_index].C_block;
   0x0000000000000029 <+41>:	mov    0x10(%rsi,%rdx,1),%r9

102	    norm = multiply_stream[stream_index].norm;
103
104	    /* Reset C(1,1) matrix accumulators */
105	    C_row[0] = _mm_setzero_ps();
   0x000000000000002e <+46>:	xorps  %xmm4,%xmm4

106	    C_row[1] = _mm_setzero_ps();
   0x0000000000000031 <+49>:	xorps  %xmm3,%xmm3

107	    C_row[2] = _mm_setzero_ps();
   0x0000000000000034 <+52>:	xorps  %xmm2,%xmm2

108	    C_row[3] = _mm_setzero_ps();
   0x0000000000000043 <+67>:	xorps  %xmm8,%xmm8

109
110	    if (norm[0]*norm[16] >= tolerance &&
   0x0000000000000037 <+55>:	movss  0x18(%rdx,%rsi,1),%xmm0
   0x000000000000003d <+61>:	movss  0x58(%rdx,%rsi,1),%xmm5
   0x0000000000000047 <+71>:	movaps %xmm0,%xmm7
   0x000000000000004a <+74>:	mulss  %xmm5,%xmm7
   0x000000000000004e <+78>:	comiss %xmm1,%xmm7
   0x0000000000000051 <+81>:	jb     0x4d2 <stream_kernel+1234>

111	        norm[1]*norm[20] >= tolerance &&
   0x0000000000000057 <+87>:	movss  0x1c(%rdx,%rsi,1),%xmm7
   0x000000000000005d <+93>:	mulss  0x68(%rdx,%rsi,1),%xmm7
   0x0000000000000063 <+99>:	comiss %xmm1,%xmm7
   0x0000000000000066 <+102>:	jb     0x4d2 <stream_kernel+1234>

112	        norm[2]*norm[24] >= tolerance &&
   0x000000000000006c <+108>:	movss  0x20(%rdx,%rsi,1),%xmm7
   0x0000000000000072 <+114>:	mulss  0x78(%rdx,%rsi,1),%xmm7
   0x0000000000000078 <+120>:	comiss %xmm1,%xmm7
   0x000000000000007b <+123>:	jb     0x4d2 <stream_kernel+1234>

113	        norm[3]*norm[28] >= tolerance)
   0x0000000000000081 <+129>:	movss  0x24(%rdx,%rsi,1),%xmm7
   0x0000000000000087 <+135>:	mulss  0x88(%rdx,%rsi,1),%xmm7
   0x0000000000000090 <+144>:	comiss %xmm1,%xmm7
   0x0000000000000093 <+147>:	jb     0x4d2 <stream_kernel+1234>

114	    {
115	      /* A(1,1)*B(1,1) = C(1,1). */
116	      for (i = 0; i < 4; i++)
117	      {
118	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
119	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
   0x0000000000000099 <+153>:	movaps (%rcx),%xmm11

120	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000000a6 <+166>:	movaps (%r8),%xmm3
   0x00000000000000b9 <+185>:	movaps 0x40(%r8),%xmm14
   0x00000000000000c8 <+200>:	movaps 0x80(%r8),%xmm9
   0x00000000000000d8 <+216>:	mulps  %xmm11,%xmm3
   0x0000000000000102 <+258>:	mulps  %xmm11,%xmm14
   0x000000000000011e <+286>:	movaps 0xc0(%r8),%xmm13
   0x0000000000000135 <+309>:	mulps  %xmm11,%xmm9
   0x0000000000000168 <+360>:	mulps  %xmm11,%xmm13

121
122	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
123	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
124	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000009d <+157>:	movaps 0x10(%rcx),%xmm10
   0x00000000000000aa <+170>:	movaps 0x10(%r8),%xmm4
   0x00000000000000be <+190>:	movaps 0x50(%r8),%xmm13
   0x00000000000000d0 <+208>:	movaps 0x90(%r8),%xmm15
   0x00000000000000dc <+220>:	mulps  %xmm10,%xmm4
   0x00000000000000e0 <+224>:	addps  %xmm3,%xmm4
   0x0000000000000106 <+262>:	mulps  %xmm10,%xmm13
   0x000000000000010a <+266>:	addps  %xmm14,%xmm13
   0x0000000000000139 <+313>:	mulps  %xmm10,%xmm15
   0x000000000000013d <+317>:	addps  %xmm9,%xmm15
   0x000000000000016c <+364>:	movaps 0xd0(%r8),%xmm11
   0x0000000000000174 <+372>:	mulps  %xmm10,%xmm11
   0x0000000000000180 <+384>:	addps  %xmm13,%xmm11

125
126	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
127	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
128	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000000a2 <+162>:	movaps 0x20(%rcx),%xmm7
   0x00000000000000af <+175>:	movaps 0x20(%r8),%xmm2
   0x00000000000000c3 <+195>:	movaps 0x60(%r8),%xmm8
   0x00000000000000eb <+235>:	mulps  %xmm7,%xmm2
   0x00000000000000ee <+238>:	addps  %xmm4,%xmm2
   0x000000000000010e <+270>:	movaps 0xa0(%r8),%xmm14
   0x0000000000000116 <+278>:	mulps  %xmm7,%xmm8
   0x000000000000011a <+282>:	addps  %xmm13,%xmm8
   0x000000000000012d <+301>:	movaps 0xe0(%r8),%xmm8
   0x0000000000000149 <+329>:	mulps  %xmm7,%xmm14
   0x000000000000014d <+333>:	addps  %xmm15,%xmm14
   0x0000000000000184 <+388>:	mulps  %xmm7,%xmm8
   0x0000000000000190 <+400>:	addps  %xmm11,%xmm8

129
130	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
131	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
132	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000000b4 <+180>:	movaps 0x30(%r8),%xmm12
   0x00000000000000e3 <+227>:	movaps 0xb0(%r8),%xmm3
   0x00000000000000f1 <+241>:	movaps 0x30(%rcx),%xmm4
   0x00000000000000f5 <+245>:	mulps  %xmm4,%xmm12
   0x00000000000000f9 <+249>:	addps  %xmm2,%xmm12
   0x00000000000000fd <+253>:	movaps 0x70(%r8),%xmm2
   0x0000000000000126 <+294>:	mulps  %xmm4,%xmm2
   0x0000000000000129 <+297>:	addps  %xmm8,%xmm2
   0x0000000000000159 <+345>:	mulps  %xmm4,%xmm3
   0x000000000000015c <+348>:	addps  %xmm14,%xmm3
   0x0000000000000178 <+376>:	movaps 0xf0(%r8),%xmm10
   0x000000000000019c <+412>:	mulps  %xmm4,%xmm10
   0x00000000000001ac <+428>:	addps  %xmm8,%xmm10

133	      }
134
135	      /* A(1,2)*B(2,1) = C(1,1). */
136	      for (i = 0; i < 4; i++)
137	      {
138	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
139	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
140	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000151 <+337>:	movaps 0x100(%rcx),%xmm15
   0x0000000000000160 <+352>:	movaps 0x180(%r8),%xmm14
   0x0000000000000188 <+392>:	movaps 0x100(%r8),%xmm7
   0x00000000000001a8 <+424>:	mulps  %xmm15,%xmm7
   0x00000000000001c8 <+456>:	addps  %xmm12,%xmm7
   0x00000000000001cc <+460>:	mulps  %xmm15,%xmm14
   0x00000000000001d0 <+464>:	movaps 0x140(%r8),%xmm12
   0x00000000000001dc <+476>:	mulps  %xmm15,%xmm12
   0x00000000000001f8 <+504>:	addps  %xmm2,%xmm12
   0x0000000000000208 <+520>:	addps  %xmm3,%xmm14
   0x000000000000024c <+588>:	movaps 0x1c0(%r8),%xmm12
   0x0000000000000264 <+612>:	mulps  %xmm15,%xmm12
   0x0000000000000296 <+662>:	addps  %xmm10,%xmm12

141
142	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
143	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
144	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000141 <+321>:	movaps 0x110(%rcx),%xmm9
   0x0000000000000194 <+404>:	movaps 0x110(%r8),%xmm11
   0x00000000000001b0 <+432>:	mulps  %xmm9,%xmm11
   0x00000000000001d8 <+472>:	addps  %xmm7,%xmm11
   0x00000000000001fc <+508>:	movaps 0x150(%r8),%xmm2
   0x000000000000020c <+524>:	movaps 0x190(%r8),%xmm3
   0x0000000000000214 <+532>:	mulps  %xmm9,%xmm2
   0x0000000000000228 <+552>:	addps  %xmm12,%xmm2
   0x0000000000000244 <+580>:	mulps  %xmm9,%xmm3
   0x0000000000000258 <+600>:	addps  %xmm14,%xmm3
   0x000000000000029a <+666>:	movaps 0x1d0(%r8),%xmm10
   0x00000000000002a2 <+674>:	mulps  %xmm9,%xmm10
   0x00000000000002ae <+686>:	addps  %xmm12,%xmm10

145
146	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
147	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
148	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000001a0 <+416>:	movaps 0x120(%r8),%xmm4
   0x00000000000001b4 <+436>:	movaps 0x120(%rcx),%xmm13
   0x00000000000001c4 <+452>:	mulps  %xmm13,%xmm4
   0x00000000000001e8 <+488>:	addps  %xmm11,%xmm4
   0x000000000000022c <+556>:	movaps 0x160(%r8),%xmm12
   0x0000000000000234 <+564>:	mulps  %xmm13,%xmm12
   0x0000000000000238 <+568>:	addps  %xmm2,%xmm12
   0x000000000000023c <+572>:	movaps 0x1a0(%r8),%xmm2
   0x0000000000000254 <+596>:	mulps  %xmm13,%xmm2
   0x0000000000000268 <+616>:	addps  %xmm3,%xmm2
   0x000000000000026b <+619>:	movaps 0x1e0(%r8),%xmm3
   0x0000000000000273 <+627>:	mulps  %xmm13,%xmm3
   0x00000000000002c6 <+710>:	addps  %xmm10,%xmm3

149
150	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
151	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
152	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000001bc <+444>:	movaps 0x130(%rcx),%xmm8
   0x00000000000001e0 <+480>:	movaps 0x1b0(%r8),%xmm7
   0x00000000000001ec <+492>:	mulps  %xmm8,%xmm7
   0x00000000000001f0 <+496>:	movaps 0x130(%r8),%xmm11
   0x0000000000000204 <+516>:	mulps  %xmm8,%xmm11
   0x0000000000000218 <+536>:	addps  %xmm4,%xmm11
   0x000000000000021c <+540>:	movaps 0x170(%r8),%xmm4
   0x0000000000000224 <+548>:	mulps  %xmm8,%xmm4
   0x0000000000000248 <+584>:	addps  %xmm12,%xmm4
   0x000000000000025c <+604>:	movaps 0x1f0(%r8),%xmm14
   0x000000000000027f <+639>:	mulps  %xmm8,%xmm14
   0x000000000000028b <+651>:	addps  %xmm2,%xmm7
   0x00000000000002d6 <+726>:	addps  %xmm3,%xmm14

153	      }
154
155	      /* A(1,3)*B(3,1) = C(1,1). */
156	      for (i = 0; i < 4; i++)
157	      {
158	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
159	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
160	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000002a6 <+678>:	movaps 0x200(%r8),%xmm9
   0x00000000000002b2 <+690>:	movaps 0x200(%rcx),%xmm12
   0x00000000000002c2 <+706>:	mulps  %xmm12,%xmm9
   0x00000000000002de <+734>:	movaps 0x240(%r8),%xmm3
   0x00000000000002e6 <+742>:	addps  %xmm11,%xmm9
   0x00000000000002ea <+746>:	mulps  %xmm12,%xmm3
   0x000000000000030a <+778>:	movaps 0x280(%r8),%xmm13
   0x0000000000000322 <+802>:	mulps  %xmm12,%xmm13
   0x0000000000000326 <+806>:	addps  %xmm4,%xmm3
   0x0000000000000335 <+821>:	addps  %xmm7,%xmm13
   0x0000000000000367 <+871>:	movaps 0x2c0(%r8),%xmm3
   0x0000000000000377 <+887>:	mulps  %xmm12,%xmm3
   0x00000000000003ba <+954>:	addps  %xmm14,%xmm3

161
162	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
163	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
164	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000277 <+631>:	movaps 0x210(%r8),%xmm13
   0x00000000000002ca <+714>:	movaps 0x210(%rcx),%xmm10
   0x00000000000002d2 <+722>:	mulps  %xmm10,%xmm13
   0x00000000000002f6 <+758>:	addps  %xmm9,%xmm13
   0x000000000000032d <+813>:	movaps 0x250(%r8),%xmm4
   0x0000000000000339 <+825>:	movaps 0x290(%r8),%xmm7
   0x0000000000000341 <+833>:	mulps  %xmm10,%xmm4
   0x0000000000000345 <+837>:	addps  %xmm3,%xmm4
   0x000000000000035f <+863>:	mulps  %xmm10,%xmm7
   0x0000000000000373 <+883>:	addps  %xmm13,%xmm7
   0x000000000000037b <+891>:	movaps 0x2d0(%r8),%xmm12
   0x0000000000000386 <+902>:	mulps  %xmm10,%xmm12
   0x00000000000003d6 <+982>:	addps  %xmm3,%xmm12

165
166	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
167	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
168	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000283 <+643>:	movaps 0x220(%r8),%xmm8
   0x00000000000002ba <+698>:	movaps 0x220(%rcx),%xmm15
   0x00000000000002da <+730>:	mulps  %xmm15,%xmm8
   0x0000000000000306 <+774>:	addps  %xmm13,%xmm8
   0x0000000000000348 <+840>:	movaps 0x260(%r8),%xmm3
   0x0000000000000350 <+848>:	mulps  %xmm15,%xmm3
   0x0000000000000354 <+852>:	addps  %xmm4,%xmm3
   0x0000000000000357 <+855>:	movaps 0x2a0(%r8),%xmm4
   0x000000000000036f <+879>:	mulps  %xmm15,%xmm4
   0x0000000000000383 <+899>:	addps  %xmm7,%xmm4
   0x000000000000039e <+926>:	movaps 0x2e0(%r8),%xmm4
   0x00000000000003a6 <+934>:	mulps  %xmm15,%xmm4
   0x00000000000003e6 <+998>:	addps  %xmm12,%xmm4

169
170	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
171	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
172	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000028e <+654>:	movaps 0x230(%r8),%xmm2
   0x00000000000002ee <+750>:	movaps 0x230(%rcx),%xmm11
   0x00000000000002fa <+762>:	mulps  %xmm11,%xmm2
   0x00000000000002fe <+766>:	movaps 0x270(%r8),%xmm9
   0x0000000000000312 <+786>:	mulps  %xmm11,%xmm9
   0x0000000000000316 <+790>:	addps  %xmm8,%xmm2
   0x000000000000031a <+794>:	movaps 0x2b0(%r8),%xmm8
   0x0000000000000329 <+809>:	mulps  %xmm11,%xmm8
   0x0000000000000363 <+867>:	addps  %xmm3,%xmm9
   0x000000000000039a <+922>:	addps  %xmm4,%xmm8
   0x00000000000003be <+958>:	movaps 0x2f0(%r8),%xmm14
   0x00000000000003c6 <+966>:	mulps  %xmm11,%xmm14
   0x00000000000003f2 <+1010>:	addps  %xmm4,%xmm14

173	      }
174
175	      /* A(1,4)*B(4,1) = C(1,1). */
176	      for (i = 0; i < 4; i++)
177	      {
178	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
179	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
180	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000038a <+906>:	movaps 0x300(%rcx),%xmm13
   0x00000000000003ca <+970>:	movaps 0x300(%r8),%xmm11
   0x00000000000003d2 <+978>:	mulps  %xmm13,%xmm11
   0x00000000000003da <+986>:	movaps 0x340(%r8),%xmm3
   0x00000000000003e2 <+994>:	mulps  %xmm13,%xmm3
   0x0000000000000402 <+1026>:	addps  %xmm2,%xmm11
   0x0000000000000412 <+1042>:	addps  %xmm9,%xmm3
   0x0000000000000416 <+1046>:	movaps 0x380(%r8),%xmm9
   0x000000000000041e <+1054>:	mulps  %xmm13,%xmm9
   0x0000000000000432 <+1074>:	addps  %xmm8,%xmm9
   0x000000000000048b <+1163>:	movaps 0x3c0(%r8),%xmm9
   0x0000000000000493 <+1171>:	mulps  %xmm13,%xmm9
   0x00000000000004c2 <+1218>:	addps  %xmm14,%xmm9

181
182	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
183	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
184	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000003ea <+1002>:	movaps 0x310(%rcx),%xmm12
   0x0000000000000406 <+1030>:	movaps 0x310(%r8),%xmm2
   0x000000000000040e <+1038>:	mulps  %xmm12,%xmm2
   0x0000000000000422 <+1058>:	addps  %xmm11,%xmm2
   0x000000000000043a <+1082>:	movaps 0x390(%r8),%xmm8
   0x0000000000000446 <+1094>:	movaps 0x350(%r8),%xmm2
   0x000000000000044e <+1102>:	mulps  %xmm12,%xmm2
   0x0000000000000456 <+1110>:	mulps  %xmm12,%xmm8
   0x000000000000045a <+1114>:	addps  %xmm3,%xmm2
   0x0000000000000487 <+1159>:	addps  %xmm9,%xmm8
   0x0000000000000497 <+1175>:	movaps 0x3d0(%r8),%xmm13
   0x00000000000004ab <+1195>:	mulps  %xmm12,%xmm13
   0x00000000000004c6 <+1222>:	addps  %xmm9,%xmm13

185
186	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
187	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
188	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000003aa <+938>:	movaps 0x320(%r8),%xmm15
   0x00000000000003b2 <+946>:	movaps 0x360(%r8),%xmm7
   0x0000000000000426 <+1062>:	movaps 0x320(%rcx),%xmm11
   0x000000000000042e <+1070>:	mulps  %xmm11,%xmm15
   0x0000000000000436 <+1078>:	mulps  %xmm11,%xmm7
   0x0000000000000442 <+1090>:	addps  %xmm2,%xmm15
   0x0000000000000469 <+1129>:	addps  %xmm2,%xmm7
   0x000000000000047b <+1147>:	movaps 0x3a0(%r8),%xmm7
   0x0000000000000483 <+1155>:	mulps  %xmm11,%xmm7
   0x000000000000049f <+1183>:	addps  %xmm8,%xmm7
   0x00000000000004af <+1199>:	movaps 0x3e0(%r8),%xmm12
   0x00000000000004b7 <+1207>:	mulps  %xmm11,%xmm12
   0x00000000000004ca <+1226>:	addps  %xmm13,%xmm12

189
190	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
191	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
192	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000392 <+914>:	movaps 0x330(%rcx),%xmm10
   0x00000000000003f6 <+1014>:	movaps 0x330(%r8),%xmm4
   0x00000000000003fe <+1022>:	mulps  %xmm10,%xmm4
   0x0000000000000452 <+1106>:	addps  %xmm15,%xmm4
   0x000000000000045d <+1117>:	movaps 0x370(%r8),%xmm3
   0x0000000000000465 <+1125>:	mulps  %xmm10,%xmm3
   0x000000000000046c <+1132>:	movaps 0x3b0(%r8),%xmm2
   0x0000000000000474 <+1140>:	mulps  %xmm10,%xmm2
   0x0000000000000478 <+1144>:	addps  %xmm7,%xmm3
   0x00000000000004a3 <+1187>:	movaps 0x3f0(%r8),%xmm8
   0x00000000000004bb <+1211>:	addps  %xmm7,%xmm2
   0x00000000000004be <+1214>:	mulps  %xmm10,%xmm8
   0x00000000000004ce <+1230>:	addps  %xmm12,%xmm8

193	      }
194	    }
195
196	    /* Store C(1,1) block. */
197	    for (i = 0; i < 4; i++)
198	    {
199	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000004d2 <+1234>:	mulps  %xmm6,%xmm4
   0x00000000000004ea <+1258>:	mulps  %xmm6,%xmm3
   0x00000000000004fa <+1274>:	mulps  %xmm6,%xmm2
   0x000000000000050a <+1290>:	mulps  %xmm6,%xmm8

200	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_11]), C_row[i]);
   0x00000000000004d5 <+1237>:	addps  (%r9),%xmm4
   0x00000000000004ed <+1261>:	addps  0x10(%r9),%xmm3
   0x00000000000004fd <+1277>:	addps  0x20(%r9),%xmm2
   0x000000000000050e <+1294>:	addps  0x30(%r9),%xmm8

201	      _mm_store_ps(&C[i*4+C_OFFSET_11], C_row[i]);
   0x00000000000004e0 <+1248>:	movaps %xmm4,(%r9)
   0x00000000000004f2 <+1266>:	movaps %xmm3,0x10(%r9)
   0x0000000000000502 <+1282>:	movaps %xmm2,0x20(%r9)
   0x0000000000000513 <+1299>:	movaps %xmm8,0x30(%r9)

202	    }
203
204	    /* Reset C(1,2) matrix accumulators */
205	    C_row[0] = _mm_setzero_ps();
   0x0000000000000518 <+1304>:	xorps  %xmm8,%xmm8

206	    C_row[1] = _mm_setzero_ps();
   0x00000000000004d9 <+1241>:	xorps  %xmm7,%xmm7

207	    C_row[2] = _mm_setzero_ps();
   0x00000000000004f7 <+1271>:	xorps  %xmm3,%xmm3

208	    C_row[3] = _mm_setzero_ps();
   0x00000000000004dc <+1244>:	xorps  %xmm10,%xmm10

209
210	    if (norm[0]*norm[17] >= tolerance &&
   0x00000000000004e4 <+1252>:	movss  0x5c(%rdx,%rsi,1),%xmm4
   0x0000000000000507 <+1287>:	movaps %xmm0,%xmm2
   0x000000000000051c <+1308>:	mulss  %xmm4,%xmm2
   0x0000000000000520 <+1312>:	comiss %xmm1,%xmm2
   0x0000000000000523 <+1315>:	jb     0x9ad <stream_kernel+2477>

211	        norm[1]*norm[21] >= tolerance &&
   0x0000000000000529 <+1321>:	movss  0x1c(%rdx,%rsi,1),%xmm2
   0x000000000000052f <+1327>:	mulss  0x6c(%rdx,%rsi,1),%xmm2
   0x0000000000000535 <+1333>:	comiss %xmm1,%xmm2
   0x0000000000000538 <+1336>:	jb     0x9ad <stream_kernel+2477>

212	        norm[2]*norm[25] >= tolerance &&
   0x000000000000053e <+1342>:	movss  0x20(%rdx,%rsi,1),%xmm2
   0x0000000000000544 <+1348>:	mulss  0x7c(%rdx,%rsi,1),%xmm2
   0x000000000000054a <+1354>:	comiss %xmm1,%xmm2
   0x000000000000054d <+1357>:	jb     0x9ad <stream_kernel+2477>

213	        norm[3]*norm[29] >= tolerance)
   0x0000000000000553 <+1363>:	movss  0x24(%rdx,%rsi,1),%xmm2
   0x0000000000000559 <+1369>:	mulss  0x8c(%rdx,%rsi,1),%xmm2
   0x0000000000000562 <+1378>:	comiss %xmm1,%xmm2
   0x0000000000000565 <+1381>:	jb     0x9ad <stream_kernel+2477>

214	    {
215	      /* A(1,1)*B(1,2) = C(1,2). */
216	      for (i = 0; i < 4; i++)
217	      {
218	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
219	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
220	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000056b <+1387>:	movaps 0x40(%rcx),%xmm11
   0x000000000000057e <+1406>:	movaps (%r8),%xmm7
   0x0000000000000591 <+1425>:	movaps 0x40(%r8),%xmm8
   0x000000000000059b <+1435>:	mulps  %xmm11,%xmm7
   0x00000000000005c9 <+1481>:	movaps 0x80(%r8),%xmm13
   0x00000000000005d1 <+1489>:	mulps  %xmm11,%xmm8
   0x00000000000005fc <+1532>:	movaps 0xc0(%r8),%xmm10
   0x0000000000000604 <+1540>:	mulps  %xmm11,%xmm13
   0x0000000000000638 <+1592>:	mulps  %xmm11,%xmm10

221
222	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
223	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
224	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000570 <+1392>:	movaps 0x50(%rcx),%xmm14
   0x0000000000000582 <+1410>:	movaps 0x10(%r8),%xmm10
   0x0000000000000596 <+1430>:	movaps 0x50(%r8),%xmm3
   0x000000000000059f <+1439>:	mulps  %xmm14,%xmm10
   0x00000000000005ab <+1451>:	addps  %xmm7,%xmm10
   0x00000000000005d5 <+1493>:	mulps  %xmm14,%xmm3
   0x00000000000005d9 <+1497>:	addps  %xmm8,%xmm3
   0x00000000000005dd <+1501>:	movaps 0x90(%r8),%xmm8
   0x0000000000000608 <+1544>:	mulps  %xmm14,%xmm8
   0x000000000000060c <+1548>:	addps  %xmm13,%xmm8
   0x000000000000063c <+1596>:	movaps 0xd0(%r8),%xmm11
   0x0000000000000644 <+1604>:	mulps  %xmm14,%xmm11
   0x0000000000000650 <+1616>:	addps  %xmm10,%xmm11

225
226	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
227	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
228	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000575 <+1397>:	movaps 0x60(%rcx),%xmm12
   0x0000000000000587 <+1415>:	movaps 0x20(%r8),%xmm13
   0x00000000000005a3 <+1443>:	movaps 0xa0(%r8),%xmm15
   0x00000000000005b4 <+1460>:	mulps  %xmm12,%xmm13
   0x00000000000005b8 <+1464>:	addps  %xmm10,%xmm13
   0x00000000000005bc <+1468>:	movaps 0x60(%r8),%xmm10
   0x00000000000005e5 <+1509>:	mulps  %xmm12,%xmm10
   0x00000000000005e9 <+1513>:	addps  %xmm3,%xmm10
   0x0000000000000618 <+1560>:	mulps  %xmm12,%xmm15
   0x000000000000061c <+1564>:	addps  %xmm8,%xmm15
   0x0000000000000648 <+1608>:	movaps 0xe0(%r8),%xmm14
   0x0000000000000654 <+1620>:	mulps  %xmm12,%xmm14
   0x0000000000000660 <+1632>:	addps  %xmm11,%xmm14

229
230	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
231	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
232	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000057a <+1402>:	movaps 0x70(%rcx),%xmm2
   0x000000000000058c <+1420>:	movaps 0x30(%r8),%xmm9
   0x00000000000005af <+1455>:	movaps 0x70(%r8),%xmm7
   0x00000000000005c1 <+1473>:	mulps  %xmm2,%xmm9
   0x00000000000005c5 <+1477>:	addps  %xmm13,%xmm9
   0x00000000000005ed <+1517>:	movaps 0xf0(%r8),%xmm3
   0x00000000000005f5 <+1525>:	mulps  %xmm2,%xmm7
   0x00000000000005f8 <+1528>:	addps  %xmm10,%xmm7
   0x0000000000000620 <+1568>:	movaps 0xb0(%r8),%xmm8
   0x0000000000000628 <+1576>:	mulps  %xmm2,%xmm8
   0x000000000000062c <+1580>:	addps  %xmm15,%xmm8
   0x000000000000066c <+1644>:	mulps  %xmm2,%xmm3
   0x000000000000067b <+1659>:	addps  %xmm14,%xmm3

233	      }
234
235	      /* A(1,2)*B(2,2) = C(1,2). */
236	      for (i = 0; i < 4; i++)
237	      {
238	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
239	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
240	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000610 <+1552>:	movaps 0x140(%rcx),%xmm13
   0x0000000000000658 <+1624>:	movaps 0x100(%r8),%xmm12
   0x0000000000000664 <+1636>:	movaps 0x140(%r8),%xmm11
   0x0000000000000677 <+1655>:	mulps  %xmm13,%xmm12
   0x000000000000067f <+1663>:	mulps  %xmm13,%xmm11
   0x0000000000000697 <+1687>:	addps  %xmm9,%xmm12
   0x00000000000006b7 <+1719>:	addps  %xmm7,%xmm11
   0x00000000000006cb <+1739>:	movaps 0x180(%r8),%xmm2
   0x00000000000006d3 <+1747>:	mulps  %xmm13,%xmm2
   0x0000000000000717 <+1815>:	addps  %xmm8,%xmm2
   0x000000000000072b <+1835>:	movaps 0x1c0(%r8),%xmm2
   0x0000000000000733 <+1843>:	mulps  %xmm13,%xmm2
   0x0000000000000773 <+1907>:	addps  %xmm3,%xmm2

241
242	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
243	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
244	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000066f <+1647>:	movaps 0x110(%r8),%xmm2
   0x0000000000000683 <+1667>:	movaps 0x150(%rcx),%xmm10
   0x0000000000000693 <+1683>:	mulps  %xmm10,%xmm2
   0x00000000000006a7 <+1703>:	addps  %xmm12,%xmm2
   0x00000000000006bb <+1723>:	movaps 0x150(%r8),%xmm7
   0x00000000000006c3 <+1731>:	mulps  %xmm10,%xmm7
   0x00000000000006e7 <+1767>:	addps  %xmm11,%xmm7
   0x000000000000071b <+1819>:	movaps 0x190(%r8),%xmm8
   0x0000000000000723 <+1827>:	mulps  %xmm10,%xmm8
   0x0000000000000727 <+1831>:	addps  %xmm2,%xmm8
   0x0000000000000737 <+1847>:	movaps 0x1d0(%r8),%xmm13
   0x000000000000073f <+1855>:	mulps  %xmm10,%xmm13
   0x000000000000078e <+1934>:	addps  %xmm2,%xmm13

245
246	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
247	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
248	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000068b <+1675>:	movaps 0x160(%rcx),%xmm14
   0x00000000000006ab <+1707>:	movaps 0x120(%r8),%xmm12
   0x00000000000006b3 <+1715>:	mulps  %xmm14,%xmm12
   0x00000000000006c7 <+1735>:	addps  %xmm2,%xmm12
   0x00000000000006db <+1755>:	movaps 0x160(%r8),%xmm12
   0x00000000000006e3 <+1763>:	mulps  %xmm14,%xmm12
   0x00000000000006f7 <+1783>:	addps  %xmm7,%xmm12
   0x00000000000006fb <+1787>:	movaps 0x1a0(%r8),%xmm7
   0x0000000000000703 <+1795>:	mulps  %xmm14,%xmm7
   0x0000000000000757 <+1879>:	addps  %xmm8,%xmm7
   0x0000000000000776 <+1910>:	movaps 0x1e0(%r8),%xmm3
   0x000000000000077e <+1918>:	mulps  %xmm14,%xmm3
   0x000000000000079d <+1949>:	addps  %xmm13,%xmm3

249
250	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
251	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
252	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000630 <+1584>:	movaps 0x170(%rcx),%xmm15
   0x000000000000069b <+1691>:	movaps 0x130(%r8),%xmm9
   0x00000000000006a3 <+1699>:	mulps  %xmm15,%xmm9
   0x00000000000006d7 <+1751>:	addps  %xmm12,%xmm9
   0x00000000000006eb <+1771>:	movaps 0x170(%r8),%xmm11
   0x00000000000006f3 <+1779>:	mulps  %xmm15,%xmm11
   0x0000000000000707 <+1799>:	addps  %xmm12,%xmm11
   0x000000000000070b <+1803>:	movaps 0x1b0(%r8),%xmm12
   0x0000000000000713 <+1811>:	mulps  %xmm15,%xmm12
   0x0000000000000743 <+1859>:	movaps 0x1f0(%r8),%xmm10
   0x000000000000074b <+1867>:	mulps  %xmm15,%xmm10
   0x0000000000000763 <+1891>:	addps  %xmm7,%xmm12
   0x00000000000007a9 <+1961>:	addps  %xmm3,%xmm10

253	      }
254
255	      /* A(1,3)*B(3,2) = C(1,2). */
256	      for (i = 0; i < 4; i++)
257	      {
258	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
259	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
260	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000075b <+1883>:	movaps 0x240(%rcx),%xmm8
   0x0000000000000767 <+1895>:	movaps 0x240(%r8),%xmm7
   0x000000000000076f <+1903>:	mulps  %xmm8,%xmm7
   0x0000000000000782 <+1922>:	movaps 0x200(%r8),%xmm14
   0x000000000000078a <+1930>:	mulps  %xmm8,%xmm14
   0x00000000000007b4 <+1972>:	addps  %xmm9,%xmm14
   0x00000000000007c4 <+1988>:	addps  %xmm11,%xmm7
   0x00000000000007d8 <+2008>:	movaps 0x280(%r8),%xmm14
   0x00000000000007e0 <+2016>:	mulps  %xmm8,%xmm14
   0x0000000000000817 <+2071>:	movaps 0x2c0(%r8),%xmm11
   0x000000000000081f <+2079>:	mulps  %xmm8,%xmm11
   0x000000000000083e <+2110>:	addps  %xmm12,%xmm14
   0x000000000000084e <+2126>:	addps  %xmm10,%xmm11

261
262	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
263	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
264	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000007ad <+1965>:	movaps 0x250(%rcx),%xmm3
   0x00000000000007b8 <+1976>:	movaps 0x210(%r8),%xmm9
   0x00000000000007c0 <+1984>:	mulps  %xmm3,%xmm9
   0x00000000000007c8 <+1992>:	movaps 0x250(%r8),%xmm11
   0x00000000000007d0 <+2000>:	mulps  %xmm3,%xmm11
   0x00000000000007d4 <+2004>:	addps  %xmm14,%xmm9
   0x00000000000007f4 <+2036>:	addps  %xmm7,%xmm11
   0x0000000000000842 <+2114>:	movaps 0x290(%r8),%xmm12
   0x000000000000084a <+2122>:	mulps  %xmm3,%xmm12
   0x0000000000000852 <+2130>:	movaps 0x2d0(%r8),%xmm10
   0x0000000000000862 <+2146>:	mulps  %xmm3,%xmm10
   0x0000000000000872 <+2162>:	addps  %xmm14,%xmm12
   0x00000000000008a2 <+2210>:	addps  %xmm11,%xmm10

265
266	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
267	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
268	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000074f <+1871>:	movaps 0x220(%r8),%xmm15
   0x0000000000000792 <+1938>:	movaps 0x260(%rcx),%xmm2
   0x0000000000000799 <+1945>:	mulps  %xmm2,%xmm15
   0x00000000000007e4 <+2020>:	addps  %xmm9,%xmm15
   0x00000000000007f8 <+2040>:	movaps 0x260(%r8),%xmm7
   0x0000000000000800 <+2048>:	mulps  %xmm2,%xmm7
   0x0000000000000813 <+2067>:	addps  %xmm11,%xmm7
   0x0000000000000823 <+2083>:	movaps 0x2e0(%r8),%xmm8
   0x000000000000082b <+2091>:	mulps  %xmm2,%xmm8
   0x0000000000000833 <+2099>:	movaps 0x2a0(%r8),%xmm7
   0x000000000000083b <+2107>:	mulps  %xmm2,%xmm7
   0x0000000000000882 <+2178>:	addps  %xmm12,%xmm7
   0x00000000000008b2 <+2226>:	addps  %xmm10,%xmm8

269
270	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
271	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
272	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000007a1 <+1953>:	movaps 0x230(%r8),%xmm13
   0x00000000000007e8 <+2024>:	movaps 0x270(%rcx),%xmm9
   0x00000000000007f0 <+2032>:	mulps  %xmm9,%xmm13
   0x0000000000000803 <+2051>:	addps  %xmm15,%xmm13
   0x0000000000000807 <+2055>:	movaps 0x270(%r8),%xmm15
   0x000000000000080f <+2063>:	mulps  %xmm9,%xmm15
   0x000000000000082f <+2095>:	addps  %xmm7,%xmm15
   0x000000000000085a <+2138>:	movaps 0x2f0(%r8),%xmm2
   0x000000000000086e <+2158>:	mulps  %xmm9,%xmm2
   0x0000000000000876 <+2166>:	movaps 0x2b0(%r8),%xmm14
   0x000000000000087e <+2174>:	mulps  %xmm9,%xmm14
   0x000000000000088e <+2190>:	addps  %xmm7,%xmm14
   0x00000000000008c2 <+2242>:	addps  %xmm8,%xmm2

273	      }
274
275	      /* A(1,4)*B(4,2) = C(1,2). */
276	      for (i = 0; i < 4; i++)
277	      {
278	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
279	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
280	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000866 <+2150>:	movaps 0x300(%r8),%xmm3
   0x000000000000089a <+2202>:	movaps 0x340(%r8),%xmm9
   0x00000000000008a6 <+2214>:	movaps 0x340(%rcx),%xmm11
   0x00000000000008ae <+2222>:	mulps  %xmm11,%xmm3
   0x00000000000008c6 <+2246>:	mulps  %xmm11,%xmm9
   0x00000000000008d2 <+2258>:	addps  %xmm13,%xmm3
   0x00000000000008e2 <+2274>:	addps  %xmm15,%xmm9
   0x0000000000000916 <+2326>:	movaps 0x380(%r8),%xmm9
   0x000000000000091e <+2334>:	mulps  %xmm11,%xmm9
   0x0000000000000936 <+2358>:	movaps 0x3c0(%r8),%xmm15
   0x000000000000093e <+2366>:	mulps  %xmm11,%xmm15
   0x0000000000000951 <+2385>:	addps  %xmm14,%xmm9
   0x0000000000000961 <+2401>:	addps  %xmm2,%xmm15

281
282	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
283	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
284	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000886 <+2182>:	movaps 0x350(%rcx),%xmm12
   0x00000000000008d6 <+2262>:	movaps 0x310(%r8),%xmm13
   0x00000000000008de <+2270>:	mulps  %xmm12,%xmm13
   0x00000000000008e6 <+2278>:	movaps 0x350(%r8),%xmm15
   0x00000000000008ee <+2286>:	mulps  %xmm12,%xmm15
   0x00000000000008f2 <+2290>:	addps  %xmm3,%xmm13
   0x0000000000000912 <+2322>:	addps  %xmm9,%xmm15
   0x0000000000000955 <+2389>:	movaps 0x390(%r8),%xmm14
   0x000000000000095d <+2397>:	mulps  %xmm12,%xmm14
   0x0000000000000965 <+2405>:	movaps 0x3d0(%r8),%xmm2
   0x000000000000096d <+2413>:	mulps  %xmm12,%xmm2
   0x0000000000000971 <+2417>:	addps  %xmm9,%xmm14
   0x0000000000000981 <+2433>:	addps  %xmm15,%xmm2

285
286	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
287	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
288	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000892 <+2194>:	movaps 0x320(%r8),%xmm7
   0x00000000000008b6 <+2230>:	movaps 0x360(%rcx),%xmm10
   0x00000000000008be <+2238>:	mulps  %xmm10,%xmm7
   0x00000000000008f6 <+2294>:	movaps 0x360(%r8),%xmm3
   0x00000000000008fe <+2302>:	mulps  %xmm10,%xmm3
   0x0000000000000902 <+2306>:	addps  %xmm13,%xmm7
   0x0000000000000932 <+2354>:	addps  %xmm15,%xmm3
   0x0000000000000975 <+2421>:	movaps 0x3a0(%r8),%xmm9
   0x000000000000097d <+2429>:	mulps  %xmm10,%xmm9
   0x0000000000000985 <+2437>:	addps  %xmm14,%xmm9
   0x000000000000098d <+2445>:	movaps 0x3e0(%r8),%xmm9
   0x0000000000000995 <+2453>:	mulps  %xmm10,%xmm9
   0x00000000000009a5 <+2469>:	addps  %xmm2,%xmm9

289
290	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
291	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
292	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000008ca <+2250>:	movaps 0x330(%r8),%xmm8
   0x0000000000000906 <+2310>:	movaps 0x370(%rcx),%xmm13
   0x000000000000090e <+2318>:	mulps  %xmm13,%xmm8
   0x0000000000000922 <+2338>:	addps  %xmm7,%xmm8
   0x0000000000000926 <+2342>:	movaps 0x370(%r8),%xmm7
   0x000000000000092e <+2350>:	mulps  %xmm13,%xmm7
   0x0000000000000942 <+2370>:	addps  %xmm3,%xmm7
   0x0000000000000945 <+2373>:	movaps 0x3b0(%r8),%xmm3
   0x000000000000094d <+2381>:	mulps  %xmm13,%xmm3
   0x0000000000000989 <+2441>:	addps  %xmm9,%xmm3
   0x0000000000000999 <+2457>:	movaps 0x3f0(%r8),%xmm10
   0x00000000000009a1 <+2465>:	mulps  %xmm13,%xmm10
   0x00000000000009a9 <+2473>:	addps  %xmm9,%xmm10

293	      }
294	    }
295
296	    /* Store C(1,2) block. */
297	    for (i = 0; i < 4; i++)
298	    {
299	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000009ad <+2477>:	mulps  %xmm6,%xmm8
   0x00000000000009c6 <+2502>:	mulps  %xmm6,%xmm7
   0x00000000000009da <+2522>:	mulps  %xmm6,%xmm3
   0x00000000000009ed <+2541>:	mulps  %xmm6,%xmm10

300	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_12]), C_row[i]);
   0x00000000000009b1 <+2481>:	addps  0x40(%r9),%xmm8
   0x00000000000009c9 <+2505>:	addps  0x50(%r9),%xmm7
   0x00000000000009dd <+2525>:	addps  0x60(%r9),%xmm3
   0x00000000000009f1 <+2545>:	addps  0x70(%r9),%xmm10

301	      _mm_store_ps(&C[i*4+C_OFFSET_12], C_row[i]);
   0x00000000000009bd <+2493>:	movaps %xmm8,0x40(%r9)
   0x00000000000009d2 <+2514>:	movaps %xmm7,0x50(%r9)
   0x00000000000009e2 <+2530>:	movaps %xmm3,0x60(%r9)
   0x00000000000009f6 <+2550>:	movaps %xmm10,0x70(%r9)

302	    }
303
304	    /* Reset C(1,3) matrix accumulators */
305	    C_row[0] = _mm_setzero_ps();
   0x00000000000009b6 <+2486>:	xorps  %xmm9,%xmm9

306	    C_row[1] = _mm_setzero_ps();
   0x00000000000009d7 <+2519>:	xorps  %xmm7,%xmm7

307	    C_row[2] = _mm_setzero_ps();
   0x00000000000009ba <+2490>:	xorps  %xmm2,%xmm2

308	    C_row[3] = _mm_setzero_ps();
   0x00000000000009c2 <+2498>:	xorps  %xmm11,%xmm11

309
310	    if (norm[0]*norm[18] >= tolerance &&
   0x00000000000009ce <+2510>:	movaps %xmm0,%xmm8
   0x00000000000009e7 <+2535>:	movss  0x60(%rdx,%rsi,1),%xmm3
   0x00000000000009fb <+2555>:	mulss  %xmm3,%xmm8
   0x0000000000000a00 <+2560>:	comiss %xmm1,%xmm8
   0x0000000000000a04 <+2564>:	jb     0xf19 <stream_kernel+3865>

311	        norm[1]*norm[22] >= tolerance &&
   0x0000000000000a0a <+2570>:	movss  0x1c(%rdx,%rsi,1),%xmm8
   0x0000000000000a11 <+2577>:	mulss  0x70(%rdx,%rsi,1),%xmm8
   0x0000000000000a18 <+2584>:	comiss %xmm1,%xmm8
   0x0000000000000a1c <+2588>:	jb     0xf19 <stream_kernel+3865>

312	        norm[2]*norm[26] >= tolerance &&
   0x0000000000000a22 <+2594>:	movss  0x20(%rdx,%rsi,1),%xmm8
   0x0000000000000a29 <+2601>:	mulss  0x80(%rdx,%rsi,1),%xmm8
   0x0000000000000a33 <+2611>:	comiss %xmm1,%xmm8
   0x0000000000000a37 <+2615>:	jb     0xf19 <stream_kernel+3865>

313	        norm[3]*norm[30] >= tolerance)
   0x0000000000000a3d <+2621>:	movss  0x24(%rdx,%rsi,1),%xmm8
   0x0000000000000a44 <+2628>:	mulss  0x90(%rdx,%rsi,1),%xmm8
   0x0000000000000a4e <+2638>:	comiss %xmm1,%xmm8
   0x0000000000000a52 <+2642>:	jb     0xf19 <stream_kernel+3865>

314	    {
315	      /* A(1,1)*B(1,3) = C(1,3). */
316	      for (i = 0; i < 4; i++)
317	      {
318	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
319	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
320	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000a58 <+2648>:	movaps (%r8),%xmm2
   0x0000000000000a6b <+2667>:	movaps 0x40(%r8),%xmm11
   0x0000000000000a7f <+2687>:	movaps 0x80(%r8),%xmm14
   0x0000000000000a8f <+2703>:	mulps  0x80(%rcx),%xmm2
   0x0000000000000aad <+2733>:	mulps  0x80(%rcx),%xmm11
   0x0000000000000acd <+2765>:	mulps  0x80(%rcx),%xmm14
   0x0000000000000b08 <+2824>:	movaps 0xc0(%r8),%xmm10
   0x0000000000000b10 <+2832>:	mulps  0x80(%rcx),%xmm10

321
322	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
323	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
324	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000a5c <+2652>:	movaps 0x10(%r8),%xmm8
   0x0000000000000a70 <+2672>:	movaps 0x50(%r8),%xmm12
   0x0000000000000a96 <+2710>:	mulps  0x90(%rcx),%xmm8
   0x0000000000000ab5 <+2741>:	mulps  0x90(%rcx),%xmm12
   0x0000000000000add <+2781>:	addps  %xmm2,%xmm8
   0x0000000000000ae1 <+2785>:	movaps 0x90(%r8),%xmm2
   0x0000000000000ae9 <+2793>:	mulps  0x90(%rcx),%xmm2
   0x0000000000000b18 <+2840>:	addps  %xmm11,%xmm12
   0x0000000000000b1c <+2844>:	movaps 0xd0(%r8),%xmm11
   0x0000000000000b24 <+2852>:	mulps  0x90(%rcx),%xmm11
   0x0000000000000b54 <+2900>:	addps  %xmm14,%xmm2
   0x0000000000000b8f <+2959>:	addps  %xmm10,%xmm11

325
326	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
327	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
328	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000a61 <+2657>:	movaps 0x20(%r8),%xmm10
   0x0000000000000a75 <+2677>:	movaps 0x60(%r8),%xmm13
   0x0000000000000a9e <+2718>:	mulps  0xa0(%rcx),%xmm10
   0x0000000000000abd <+2749>:	mulps  0xa0(%rcx),%xmm13
   0x0000000000000af0 <+2800>:	addps  %xmm8,%xmm10
   0x0000000000000af4 <+2804>:	movaps 0xa0(%r8),%xmm8
   0x0000000000000afc <+2812>:	mulps  0xa0(%rcx),%xmm8
   0x0000000000000b2c <+2860>:	addps  %xmm12,%xmm13
   0x0000000000000b30 <+2864>:	movaps 0xe0(%r8),%xmm12
   0x0000000000000b38 <+2872>:	mulps  0xa0(%rcx),%xmm12
   0x0000000000000b68 <+2920>:	addps  %xmm2,%xmm8
   0x0000000000000ba3 <+2979>:	addps  %xmm11,%xmm12

329
330	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
331	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
332	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000a66 <+2662>:	movaps 0x30(%r8),%xmm7
   0x0000000000000a7a <+2682>:	movaps 0x70(%r8),%xmm9
   0x0000000000000a87 <+2695>:	movaps 0xb0(%r8),%xmm15
   0x0000000000000aa6 <+2726>:	mulps  0xb0(%rcx),%xmm7
   0x0000000000000ac5 <+2757>:	mulps  0xb0(%rcx),%xmm9
   0x0000000000000ad5 <+2773>:	mulps  0xb0(%rcx),%xmm15
   0x0000000000000b04 <+2820>:	addps  %xmm10,%xmm7
   0x0000000000000b40 <+2880>:	addps  %xmm13,%xmm9
   0x0000000000000b58 <+2904>:	movaps 0xf0(%r8),%xmm14
   0x0000000000000b60 <+2912>:	mulps  0xb0(%rcx),%xmm14
   0x0000000000000b7b <+2939>:	addps  %xmm8,%xmm15
   0x0000000000000bb7 <+2999>:	addps  %xmm12,%xmm14

333	      }
334
335	      /* A(1,2)*B(2,3) = C(1,3). */
336	      for (i = 0; i < 4; i++)
337	      {
338	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
339	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
340	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000b6c <+2924>:	movaps 0x100(%r8),%xmm2
   0x0000000000000b74 <+2932>:	mulps  0x180(%rcx),%xmm2
   0x0000000000000b7f <+2943>:	movaps 0x140(%r8),%xmm8
   0x0000000000000b87 <+2951>:	mulps  0x180(%rcx),%xmm8
   0x0000000000000bcb <+3019>:	addps  %xmm7,%xmm2
   0x0000000000000bdd <+3037>:	addps  %xmm9,%xmm8
   0x0000000000000c2d <+3117>:	movaps 0x180(%r8),%xmm2
   0x0000000000000c35 <+3125>:	mulps  0x180(%rcx),%xmm2
   0x0000000000000c54 <+3156>:	movaps 0x1c0(%r8),%xmm7
   0x0000000000000c5c <+3164>:	mulps  0x180(%rcx),%xmm7
   0x0000000000000c63 <+3171>:	addps  %xmm15,%xmm2
   0x0000000000000c77 <+3191>:	addps  %xmm14,%xmm7

341
342	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
343	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
344	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000bce <+3022>:	movaps 0x110(%r8),%xmm7
   0x0000000000000bd6 <+3030>:	mulps  0x190(%rcx),%xmm7
   0x0000000000000be1 <+3041>:	movaps 0x150(%r8),%xmm9
   0x0000000000000be9 <+3049>:	mulps  0x190(%rcx),%xmm9
   0x0000000000000bf1 <+3057>:	addps  %xmm2,%xmm7
   0x0000000000000c03 <+3075>:	addps  %xmm8,%xmm9
   0x0000000000000c67 <+3175>:	movaps 0x190(%r8),%xmm15
   0x0000000000000c6f <+3183>:	mulps  0x190(%rcx),%xmm15
   0x0000000000000c8b <+3211>:	addps  %xmm2,%xmm15
   0x0000000000000c8f <+3215>:	movaps 0x1d0(%r8),%xmm2
   0x0000000000000c97 <+3223>:	mulps  0x190(%rcx),%xmm2
   0x0000000000000cc6 <+3270>:	addps  %xmm7,%xmm2

345
346	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
347	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
348	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000bf4 <+3060>:	movaps 0x120(%r8),%xmm2
   0x0000000000000bfc <+3068>:	mulps  0x1a0(%rcx),%xmm2
   0x0000000000000c07 <+3079>:	movaps 0x1e0(%r8),%xmm8
   0x0000000000000c0f <+3087>:	mulps  0x1a0(%rcx),%xmm8
   0x0000000000000c17 <+3095>:	addps  %xmm7,%xmm2
   0x0000000000000c1a <+3098>:	movaps 0x160(%r8),%xmm7
   0x0000000000000c22 <+3106>:	mulps  0x1a0(%rcx),%xmm7
   0x0000000000000c3c <+3132>:	addps  %xmm9,%xmm7
   0x0000000000000c40 <+3136>:	movaps 0x1a0(%r8),%xmm9
   0x0000000000000c48 <+3144>:	mulps  0x1a0(%rcx),%xmm9
   0x0000000000000c9e <+3230>:	addps  %xmm15,%xmm9
   0x0000000000000cd8 <+3288>:	addps  %xmm2,%xmm8

349
350	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
351	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
352	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000b44 <+2884>:	movaps 0x130(%r8),%xmm13
   0x0000000000000b4c <+2892>:	mulps  0x1b0(%rcx),%xmm13
   0x0000000000000b93 <+2963>:	movaps 0x1f0(%r8),%xmm10
   0x0000000000000b9b <+2971>:	mulps  0x1b0(%rcx),%xmm10
   0x0000000000000ba7 <+2983>:	movaps 0x1b0(%r8),%xmm11
   0x0000000000000baf <+2991>:	mulps  0x1b0(%rcx),%xmm11
   0x0000000000000bbb <+3003>:	movaps 0x170(%r8),%xmm12
   0x0000000000000bc3 <+3011>:	mulps  0x1b0(%rcx),%xmm12
   0x0000000000000c29 <+3113>:	addps  %xmm2,%xmm13
   0x0000000000000c50 <+3152>:	addps  %xmm7,%xmm12
   0x0000000000000cb2 <+3250>:	addps  %xmm9,%xmm11
   0x0000000000000ceb <+3307>:	addps  %xmm8,%xmm10

353	      }
354
355	      /* A(1,3)*B(3,3) = C(1,3). */
356	      for (i = 0; i < 4; i++)
357	      {
358	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
359	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
360	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000c7b <+3195>:	movaps 0x200(%r8),%xmm14
   0x0000000000000c83 <+3203>:	mulps  0x280(%rcx),%xmm14
   0x0000000000000cdc <+3292>:	movaps 0x240(%r8),%xmm2
   0x0000000000000ce4 <+3300>:	mulps  0x280(%rcx),%xmm2
   0x0000000000000cff <+3327>:	addps  %xmm13,%xmm14
   0x0000000000000d13 <+3347>:	addps  %xmm12,%xmm2
   0x0000000000000d2b <+3371>:	movaps 0x280(%r8),%xmm14
   0x0000000000000d33 <+3379>:	mulps  0x280(%rcx),%xmm14
   0x0000000000000d79 <+3449>:	movaps 0x2c0(%r8),%xmm8
   0x0000000000000d81 <+3457>:	mulps  0x280(%rcx),%xmm8
   0x0000000000000d9d <+3485>:	addps  %xmm11,%xmm14
   0x0000000000000db1 <+3505>:	addps  %xmm10,%xmm8

361
362	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
363	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
364	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000cef <+3311>:	movaps 0x250(%r8),%xmm8
   0x0000000000000cf7 <+3319>:	mulps  0x290(%rcx),%xmm8
   0x0000000000000d03 <+3331>:	movaps 0x210(%r8),%xmm13
   0x0000000000000d0b <+3339>:	mulps  0x290(%rcx),%xmm13
   0x0000000000000d27 <+3367>:	addps  %xmm14,%xmm13
   0x0000000000000d62 <+3426>:	addps  %xmm2,%xmm8
   0x0000000000000da1 <+3489>:	movaps 0x290(%r8),%xmm11
   0x0000000000000da9 <+3497>:	mulps  0x290(%rcx),%xmm11
   0x0000000000000db5 <+3509>:	movaps 0x2d0(%r8),%xmm10
   0x0000000000000dbd <+3517>:	mulps  0x290(%rcx),%xmm10
   0x0000000000000dc5 <+3525>:	addps  %xmm14,%xmm11
   0x0000000000000e01 <+3585>:	addps  %xmm8,%xmm10

365
366	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
367	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
368	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000ca2 <+3234>:	movaps 0x2a0(%r8),%xmm15
   0x0000000000000caa <+3242>:	mulps  0x2a0(%rcx),%xmm15
   0x0000000000000cc9 <+3273>:	movaps 0x220(%r8),%xmm7
   0x0000000000000cd1 <+3281>:	mulps  0x2a0(%rcx),%xmm7
   0x0000000000000d17 <+3351>:	movaps 0x260(%r8),%xmm12
   0x0000000000000d1f <+3359>:	mulps  0x2a0(%rcx),%xmm12
   0x0000000000000d3b <+3387>:	addps  %xmm13,%xmm7
   0x0000000000000d75 <+3445>:	addps  %xmm8,%xmm12
   0x0000000000000dd9 <+3545>:	addps  %xmm11,%xmm15
   0x0000000000000ddd <+3549>:	movaps 0x2e0(%r8),%xmm11
   0x0000000000000de5 <+3557>:	mulps  0x2a0(%rcx),%xmm11
   0x0000000000000e15 <+3605>:	addps  %xmm10,%xmm11

369
370	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
371	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
372	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000cb6 <+3254>:	movaps 0x230(%r8),%xmm9
   0x0000000000000cbe <+3262>:	mulps  0x2b0(%rcx),%xmm9
   0x0000000000000d4f <+3407>:	addps  %xmm7,%xmm9
   0x0000000000000d53 <+3411>:	movaps 0x270(%r8),%xmm7
   0x0000000000000d5b <+3419>:	mulps  0x2b0(%rcx),%xmm7
   0x0000000000000d66 <+3430>:	movaps 0x2b0(%r8),%xmm2
   0x0000000000000d6e <+3438>:	mulps  0x2b0(%rcx),%xmm2
   0x0000000000000d89 <+3465>:	addps  %xmm12,%xmm7
   0x0000000000000ded <+3565>:	addps  %xmm15,%xmm2
   0x0000000000000e05 <+3589>:	movaps 0x2f0(%r8),%xmm8
   0x0000000000000e0d <+3597>:	mulps  0x2b0(%rcx),%xmm8
   0x0000000000000e29 <+3625>:	addps  %xmm11,%xmm8

373	      }
374
375	      /* A(1,4)*B(4,3) = C(1,3). */
376	      for (i = 0; i < 4; i++)
377	      {
378	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
379	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
380	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d8d <+3469>:	movaps 0x300(%r8),%xmm12
   0x0000000000000d95 <+3477>:	mulps  0x380(%rcx),%xmm12
   0x0000000000000dc9 <+3529>:	movaps 0x340(%r8),%xmm14
   0x0000000000000dd1 <+3537>:	mulps  0x380(%rcx),%xmm14
   0x0000000000000e19 <+3609>:	movaps 0x380(%r8),%xmm10
   0x0000000000000e21 <+3617>:	mulps  0x380(%rcx),%xmm10
   0x0000000000000e3d <+3645>:	addps  %xmm9,%xmm12
   0x0000000000000e51 <+3665>:	addps  %xmm7,%xmm14
   0x0000000000000e68 <+3688>:	movaps 0x3c0(%r8),%xmm12
   0x0000000000000e70 <+3696>:	mulps  0x380(%rcx),%xmm12
   0x0000000000000ea7 <+3751>:	addps  %xmm2,%xmm10
   0x0000000000000ebe <+3774>:	addps  %xmm8,%xmm12

381
382	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
383	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
384	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000e41 <+3649>:	movaps 0x310(%r8),%xmm9
   0x0000000000000e49 <+3657>:	mulps  0x390(%rcx),%xmm9
   0x0000000000000e55 <+3669>:	movaps 0x350(%r8),%xmm7
   0x0000000000000e5d <+3677>:	mulps  0x390(%rcx),%xmm7
   0x0000000000000e64 <+3684>:	addps  %xmm12,%xmm9
   0x0000000000000e8c <+3724>:	addps  %xmm14,%xmm7
   0x0000000000000eab <+3755>:	movaps 0x390(%r8),%xmm2
   0x0000000000000eb3 <+3763>:	mulps  0x390(%rcx),%xmm2
   0x0000000000000ec2 <+3778>:	movaps 0x3d0(%r8),%xmm8
   0x0000000000000eca <+3786>:	mulps  0x390(%rcx),%xmm8
   0x0000000000000ed2 <+3794>:	addps  %xmm10,%xmm2
   0x0000000000000ef9 <+3833>:	addps  %xmm12,%xmm8

385
386	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
387	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
388	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000d3f <+3391>:	movaps 0x320(%r8),%xmm13
   0x0000000000000d47 <+3399>:	mulps  0x3a0(%rcx),%xmm13
   0x0000000000000df1 <+3569>:	movaps 0x360(%r8),%xmm15
   0x0000000000000df9 <+3577>:	mulps  0x3a0(%rcx),%xmm15
   0x0000000000000e2d <+3629>:	movaps 0x3a0(%r8),%xmm11
   0x0000000000000e35 <+3637>:	mulps  0x3a0(%rcx),%xmm11
   0x0000000000000e78 <+3704>:	addps  %xmm9,%xmm13
   0x0000000000000e90 <+3728>:	addps  %xmm7,%xmm15
   0x0000000000000ed6 <+3798>:	movaps 0x3e0(%r8),%xmm10
   0x0000000000000ede <+3806>:	mulps  0x3a0(%rcx),%xmm10
   0x0000000000000ee6 <+3814>:	addps  %xmm2,%xmm11
   0x0000000000000efd <+3837>:	addps  %xmm8,%xmm10

389
390	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
391	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
392	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000e7c <+3708>:	movaps 0x330(%r8),%xmm9
   0x0000000000000e84 <+3716>:	mulps  0x3b0(%rcx),%xmm9
   0x0000000000000e94 <+3732>:	movaps 0x370(%r8),%xmm7
   0x0000000000000e9c <+3740>:	mulps  0x3b0(%rcx),%xmm7
   0x0000000000000ea3 <+3747>:	addps  %xmm13,%xmm9
   0x0000000000000eba <+3770>:	addps  %xmm15,%xmm7
   0x0000000000000eea <+3818>:	movaps 0x3b0(%r8),%xmm2
   0x0000000000000ef2 <+3826>:	mulps  0x3b0(%rcx),%xmm2
   0x0000000000000f01 <+3841>:	addps  %xmm11,%xmm2
   0x0000000000000f05 <+3845>:	movaps 0x3f0(%r8),%xmm11
   0x0000000000000f0d <+3853>:	mulps  0x3b0(%rcx),%xmm11
   0x0000000000000f15 <+3861>:	addps  %xmm10,%xmm11

393	      }
394	    }
395
396	    /* Store C(1,3) block. */
397	    for (i = 0; i < 4; i++)
398	    {
399	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000000f19 <+3865>:	mulps  %xmm6,%xmm9
   0x0000000000000f35 <+3893>:	mulps  %xmm6,%xmm7
   0x0000000000000f4b <+3915>:	mulps  %xmm6,%xmm2
   0x0000000000000f64 <+3940>:	mulps  %xmm6,%xmm11

400	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_13]), C_row[i]);
   0x0000000000000f1d <+3869>:	addps  0x80(%r9),%xmm9
   0x0000000000000f38 <+3896>:	addps  0x90(%r9),%xmm7
   0x0000000000000f4e <+3918>:	addps  0xa0(%r9),%xmm2
   0x0000000000000f68 <+3944>:	addps  0xb0(%r9),%xmm11

401	      _mm_store_ps(&C[i*4+C_OFFSET_13], C_row[i]);
   0x0000000000000f29 <+3881>:	movaps %xmm9,0x80(%r9)
   0x0000000000000f40 <+3904>:	movaps %xmm7,0x90(%r9)
   0x0000000000000f56 <+3926>:	movaps %xmm2,0xa0(%r9)
   0x0000000000000f70 <+3952>:	movaps %xmm11,0xb0(%r9)

402	    }
403
404	    /* Reset C(1,4) matrix accumulators */
405	    C_row[0] = _mm_setzero_ps();
   0x0000000000000f31 <+3889>:	xorps  %xmm9,%xmm9

406	    C_row[1] = _mm_setzero_ps();
   0x0000000000000f25 <+3877>:	xorps  %xmm8,%xmm8

407	    C_row[2] = _mm_setzero_ps();
   0x0000000000000f48 <+3912>:	xorps  %xmm7,%xmm7

408	    C_row[3] = _mm_setzero_ps();
   0x0000000000000f78 <+3960>:	xorps  %xmm11,%xmm11

409
410	    if (norm[0]*norm[19] >= tolerance &&
   0x0000000000000f5e <+3934>:	movss  0x64(%rdx,%rsi,1),%xmm2
   0x0000000000000f7c <+3964>:	mulss  %xmm2,%xmm0
   0x0000000000000f80 <+3968>:	comiss %xmm1,%xmm0
   0x0000000000000f83 <+3971>:	jb     0x148f <stream_kernel+5263>

411	        norm[1]*norm[23] >= tolerance &&
   0x0000000000000f89 <+3977>:	movss  0x1c(%rdx,%rsi,1),%xmm0
   0x0000000000000f8f <+3983>:	mulss  0x74(%rdx,%rsi,1),%xmm0
   0x0000000000000f95 <+3989>:	comiss %xmm1,%xmm0
   0x0000000000000f98 <+3992>:	jb     0x148f <stream_kernel+5263>

412	        norm[2]*norm[27] >= tolerance &&
   0x0000000000000f9e <+3998>:	movss  0x20(%rdx,%rsi,1),%xmm0
   0x0000000000000fa4 <+4004>:	mulss  0x84(%rdx,%rsi,1),%xmm0
   0x0000000000000fad <+4013>:	comiss %xmm1,%xmm0
   0x0000000000000fb0 <+4016>:	jb     0x148f <stream_kernel+5263>

413	        norm[3]*norm[31] >= tolerance)
   0x0000000000000fb6 <+4022>:	movss  0x24(%rdx,%rsi,1),%xmm0
   0x0000000000000fbc <+4028>:	mulss  0x94(%rdx,%rsi,1),%xmm0
   0x0000000000000fc5 <+4037>:	comiss %xmm1,%xmm0
   0x0000000000000fc8 <+4040>:	jb     0x148f <stream_kernel+5263>

414	    {
415	      /* A(1,1)*B(1,4) = C(1,4). */
416	      for (i = 0; i < 4; i++)
417	      {
418	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
419	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
420	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000fce <+4046>:	movaps (%r8),%xmm7
   0x0000000000000fe1 <+4065>:	movaps 0x40(%r8),%xmm11
   0x0000000000000ff5 <+4085>:	movaps 0x80(%r8),%xmm14
   0x0000000000001005 <+4101>:	mulps  0xc0(%rcx),%xmm7
   0x0000000000001023 <+4131>:	mulps  0xc0(%rcx),%xmm11
   0x0000000000001043 <+4163>:	mulps  0xc0(%rcx),%xmm14
   0x000000000000107c <+4220>:	movaps 0xc0(%r8),%xmm10
   0x0000000000001084 <+4228>:	mulps  0xc0(%rcx),%xmm10

421
422	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
423	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
424	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000fd2 <+4050>:	movaps 0x10(%r8),%xmm0
   0x0000000000000fe6 <+4070>:	movaps 0x50(%r8),%xmm12
   0x000000000000100c <+4108>:	mulps  0xd0(%rcx),%xmm0
   0x000000000000102b <+4139>:	mulps  0xd0(%rcx),%xmm12
   0x0000000000001053 <+4179>:	addps  %xmm7,%xmm0
   0x0000000000001056 <+4182>:	movaps 0x90(%r8),%xmm7
   0x000000000000105e <+4190>:	mulps  0xd0(%rcx),%xmm7
   0x000000000000108c <+4236>:	addps  %xmm11,%xmm12
   0x0000000000001090 <+4240>:	movaps 0xd0(%r8),%xmm11
   0x0000000000001098 <+4248>:	mulps  0xd0(%rcx),%xmm11
   0x00000000000010c8 <+4296>:	addps  %xmm14,%xmm7
   0x0000000000001101 <+4353>:	addps  %xmm10,%xmm11

425
426	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
427	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
428	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000fd7 <+4055>:	movaps 0x20(%r8),%xmm10
   0x0000000000000feb <+4075>:	movaps 0x60(%r8),%xmm13
   0x0000000000001013 <+4115>:	mulps  0xe0(%rcx),%xmm10
   0x0000000000001033 <+4147>:	mulps  0xe0(%rcx),%xmm13
   0x0000000000001065 <+4197>:	addps  %xmm0,%xmm10
   0x0000000000001069 <+4201>:	movaps 0xa0(%r8),%xmm0
   0x0000000000001071 <+4209>:	mulps  0xe0(%rcx),%xmm0
   0x00000000000010a0 <+4256>:	addps  %xmm12,%xmm13
   0x00000000000010a4 <+4260>:	movaps 0xe0(%r8),%xmm12
   0x00000000000010ac <+4268>:	mulps  0xe0(%rcx),%xmm12
   0x00000000000010dc <+4316>:	addps  %xmm7,%xmm0
   0x0000000000001115 <+4373>:	addps  %xmm11,%xmm12

429
430	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
431	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
432	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000000fdc <+4060>:	movaps 0x30(%r8),%xmm8
   0x0000000000000ff0 <+4080>:	movaps 0x70(%r8),%xmm9
   0x0000000000000ffd <+4093>:	movaps 0xb0(%r8),%xmm15
   0x000000000000101b <+4123>:	mulps  0xf0(%rcx),%xmm8
   0x000000000000103b <+4155>:	mulps  0xf0(%rcx),%xmm9
   0x000000000000104b <+4171>:	mulps  0xf0(%rcx),%xmm15
   0x0000000000001078 <+4216>:	addps  %xmm10,%xmm8
   0x00000000000010b4 <+4276>:	addps  %xmm13,%xmm9
   0x00000000000010cc <+4300>:	movaps 0xf0(%r8),%xmm14
   0x00000000000010d4 <+4308>:	mulps  0xf0(%rcx),%xmm14
   0x00000000000010ee <+4334>:	addps  %xmm0,%xmm15
   0x0000000000001129 <+4393>:	addps  %xmm12,%xmm14

433	      }
434
435	      /* A(1,2)*B(2,4) = C(1,4). */
436	      for (i = 0; i < 4; i++)
437	      {
438	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
439	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
440	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000010df <+4319>:	movaps 0x100(%r8),%xmm7
   0x00000000000010e7 <+4327>:	mulps  0x1c0(%rcx),%xmm7
   0x00000000000010f2 <+4338>:	movaps 0x140(%r8),%xmm0
   0x00000000000010fa <+4346>:	mulps  0x1c0(%rcx),%xmm0
   0x000000000000113d <+4413>:	addps  %xmm8,%xmm7
   0x0000000000001151 <+4433>:	addps  %xmm9,%xmm0
   0x00000000000011a3 <+4515>:	movaps 0x180(%r8),%xmm7
   0x00000000000011ab <+4523>:	mulps  0x1c0(%rcx),%xmm7
   0x00000000000011ca <+4554>:	movaps 0x1c0(%r8),%xmm8
   0x00000000000011d2 <+4562>:	mulps  0x1c0(%rcx),%xmm8
   0x00000000000011da <+4570>:	addps  %xmm15,%xmm7
   0x00000000000011ee <+4590>:	addps  %xmm14,%xmm8

441
442	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
443	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
444	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001141 <+4417>:	movaps 0x110(%r8),%xmm8
   0x0000000000001149 <+4425>:	mulps  0x1d0(%rcx),%xmm8
   0x0000000000001155 <+4437>:	movaps 0x150(%r8),%xmm9
   0x000000000000115d <+4445>:	mulps  0x1d0(%rcx),%xmm9
   0x0000000000001165 <+4453>:	addps  %xmm7,%xmm8
   0x0000000000001178 <+4472>:	addps  %xmm0,%xmm9
   0x00000000000011de <+4574>:	movaps 0x190(%r8),%xmm15
   0x00000000000011e6 <+4582>:	mulps  0x1d0(%rcx),%xmm15
   0x0000000000001202 <+4610>:	addps  %xmm7,%xmm15
   0x0000000000001206 <+4614>:	movaps 0x1d0(%r8),%xmm7
   0x000000000000120e <+4622>:	mulps  0x1d0(%rcx),%xmm7
   0x000000000000123d <+4669>:	addps  %xmm8,%xmm7

445
446	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
447	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
448	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001169 <+4457>:	movaps 0x120(%r8),%xmm7
   0x0000000000001171 <+4465>:	mulps  0x1e0(%rcx),%xmm7
   0x000000000000117c <+4476>:	movaps 0x1e0(%r8),%xmm0
   0x0000000000001184 <+4484>:	mulps  0x1e0(%rcx),%xmm0
   0x000000000000118b <+4491>:	addps  %xmm8,%xmm7
   0x000000000000118f <+4495>:	movaps 0x160(%r8),%xmm8
   0x0000000000001197 <+4503>:	mulps  0x1e0(%rcx),%xmm8
   0x00000000000011b2 <+4530>:	addps  %xmm9,%xmm8
   0x00000000000011b6 <+4534>:	movaps 0x1a0(%r8),%xmm9
   0x00000000000011be <+4542>:	mulps  0x1e0(%rcx),%xmm9
   0x0000000000001215 <+4629>:	addps  %xmm15,%xmm9
   0x0000000000001251 <+4689>:	addps  %xmm7,%xmm0

449
450	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
451	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
452	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000010b8 <+4280>:	movaps 0x130(%r8),%xmm13
   0x00000000000010c0 <+4288>:	mulps  0x1f0(%rcx),%xmm13
   0x0000000000001105 <+4357>:	movaps 0x1f0(%r8),%xmm10
   0x000000000000110d <+4365>:	mulps  0x1f0(%rcx),%xmm10
   0x0000000000001119 <+4377>:	movaps 0x1b0(%r8),%xmm11
   0x0000000000001121 <+4385>:	mulps  0x1f0(%rcx),%xmm11
   0x000000000000112d <+4397>:	movaps 0x170(%r8),%xmm12
   0x0000000000001135 <+4405>:	mulps  0x1f0(%rcx),%xmm12
   0x000000000000119f <+4511>:	addps  %xmm7,%xmm13
   0x00000000000011c6 <+4550>:	addps  %xmm8,%xmm12
   0x0000000000001229 <+4649>:	addps  %xmm9,%xmm11
   0x0000000000001263 <+4707>:	addps  %xmm0,%xmm10

453	      }
454
455	      /* A(1,3)*B(3,4) = C(1,4). */
456	      for (i = 0; i < 4; i++)
457	      {
458	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_13]);
459	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
460	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000011f2 <+4594>:	movaps 0x200(%r8),%xmm14
   0x00000000000011fa <+4602>:	mulps  0x2c0(%rcx),%xmm14
   0x0000000000001254 <+4692>:	movaps 0x240(%r8),%xmm7
   0x000000000000125c <+4700>:	mulps  0x2c0(%rcx),%xmm7
   0x0000000000001276 <+4726>:	addps  %xmm13,%xmm14
   0x000000000000128a <+4746>:	addps  %xmm12,%xmm7
   0x00000000000012a2 <+4770>:	movaps 0x280(%r8),%xmm14
   0x00000000000012aa <+4778>:	mulps  0x2c0(%rcx),%xmm14
   0x00000000000012f0 <+4848>:	movaps 0x2c0(%r8),%xmm0
   0x00000000000012f8 <+4856>:	mulps  0x2c0(%rcx),%xmm0
   0x0000000000001313 <+4883>:	addps  %xmm11,%xmm14
   0x0000000000001327 <+4903>:	addps  %xmm10,%xmm0

461
462	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_13]);
463	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
464	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001267 <+4711>:	movaps 0x250(%r8),%xmm0
   0x000000000000126f <+4719>:	mulps  0x2d0(%rcx),%xmm0
   0x000000000000127a <+4730>:	movaps 0x210(%r8),%xmm13
   0x0000000000001282 <+4738>:	mulps  0x2d0(%rcx),%xmm13
   0x000000000000129e <+4766>:	addps  %xmm14,%xmm13
   0x00000000000012da <+4826>:	addps  %xmm7,%xmm0
   0x0000000000001317 <+4887>:	movaps 0x290(%r8),%xmm11
   0x000000000000131f <+4895>:	mulps  0x2d0(%rcx),%xmm11
   0x000000000000132b <+4907>:	movaps 0x2d0(%r8),%xmm10
   0x0000000000001333 <+4915>:	mulps  0x2d0(%rcx),%xmm10
   0x000000000000133b <+4923>:	addps  %xmm14,%xmm11
   0x0000000000001377 <+4983>:	addps  %xmm0,%xmm10

465
466	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_13]);
467	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
468	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001219 <+4633>:	movaps 0x2a0(%r8),%xmm15
   0x0000000000001221 <+4641>:	mulps  0x2e0(%rcx),%xmm15
   0x0000000000001241 <+4673>:	movaps 0x220(%r8),%xmm8
   0x0000000000001249 <+4681>:	mulps  0x2e0(%rcx),%xmm8
   0x000000000000128e <+4750>:	movaps 0x260(%r8),%xmm12
   0x0000000000001296 <+4758>:	mulps  0x2e0(%rcx),%xmm12
   0x00000000000012b2 <+4786>:	addps  %xmm13,%xmm8
   0x00000000000012ec <+4844>:	addps  %xmm0,%xmm12
   0x000000000000134f <+4943>:	addps  %xmm11,%xmm15
   0x0000000000001353 <+4947>:	movaps 0x2e0(%r8),%xmm11
   0x000000000000135b <+4955>:	mulps  0x2e0(%rcx),%xmm11
   0x000000000000138a <+5002>:	addps  %xmm10,%xmm11

469
470	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_13]);
471	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
472	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000122d <+4653>:	movaps 0x230(%r8),%xmm9
   0x0000000000001235 <+4661>:	mulps  0x2f0(%rcx),%xmm9
   0x00000000000012c6 <+4806>:	addps  %xmm8,%xmm9
   0x00000000000012ca <+4810>:	movaps 0x270(%r8),%xmm8
   0x00000000000012d2 <+4818>:	mulps  0x2f0(%rcx),%xmm8
   0x00000000000012dd <+4829>:	movaps 0x2b0(%r8),%xmm7
   0x00000000000012e5 <+4837>:	mulps  0x2f0(%rcx),%xmm7
   0x00000000000012ff <+4863>:	addps  %xmm12,%xmm8
   0x0000000000001363 <+4963>:	addps  %xmm15,%xmm7
   0x000000000000137b <+4987>:	movaps 0x2f0(%r8),%xmm0
   0x0000000000001383 <+4995>:	mulps  0x2f0(%rcx),%xmm0
   0x000000000000139e <+5022>:	addps  %xmm11,%xmm0

473	      }
474
475	      /* A(1,4)*B(4,4) = C(1,4). */
476	      for (i = 0; i < 4; i++)
477	      {
478	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_14]);
479	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
480	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001303 <+4867>:	movaps 0x300(%r8),%xmm12
   0x000000000000130b <+4875>:	mulps  0x3c0(%rcx),%xmm12
   0x000000000000133f <+4927>:	movaps 0x340(%r8),%xmm14
   0x0000000000001347 <+4935>:	mulps  0x3c0(%rcx),%xmm14
   0x000000000000138e <+5006>:	movaps 0x380(%r8),%xmm10
   0x0000000000001396 <+5014>:	mulps  0x3c0(%rcx),%xmm10
   0x00000000000013b2 <+5042>:	addps  %xmm9,%xmm12
   0x00000000000013c6 <+5062>:	addps  %xmm8,%xmm14
   0x00000000000013de <+5086>:	movaps 0x3c0(%r8),%xmm12
   0x00000000000013e6 <+5094>:	mulps  0x3c0(%rcx),%xmm12
   0x000000000000141e <+5150>:	addps  %xmm7,%xmm10
   0x0000000000001435 <+5173>:	addps  %xmm0,%xmm12

481
482	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_14]);
483	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
484	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000013b6 <+5046>:	movaps 0x310(%r8),%xmm9
   0x00000000000013be <+5054>:	mulps  0x3d0(%rcx),%xmm9
   0x00000000000013ca <+5066>:	movaps 0x350(%r8),%xmm8
   0x00000000000013d2 <+5074>:	mulps  0x3d0(%rcx),%xmm8
   0x00000000000013da <+5082>:	addps  %xmm12,%xmm9
   0x0000000000001402 <+5122>:	addps  %xmm14,%xmm8
   0x0000000000001422 <+5154>:	movaps 0x390(%r8),%xmm7
   0x000000000000142a <+5162>:	mulps  0x3d0(%rcx),%xmm7
   0x0000000000001439 <+5177>:	movaps 0x3d0(%r8),%xmm0
   0x0000000000001441 <+5185>:	mulps  0x3d0(%rcx),%xmm0
   0x0000000000001448 <+5192>:	addps  %xmm10,%xmm7
   0x000000000000146f <+5231>:	addps  %xmm12,%xmm0

485
486	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_14]);
487	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
488	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000012b6 <+4790>:	movaps 0x320(%r8),%xmm13
   0x00000000000012be <+4798>:	mulps  0x3e0(%rcx),%xmm13
   0x0000000000001367 <+4967>:	movaps 0x360(%r8),%xmm15
   0x000000000000136f <+4975>:	mulps  0x3e0(%rcx),%xmm15
   0x00000000000013a2 <+5026>:	movaps 0x3a0(%r8),%xmm11
   0x00000000000013aa <+5034>:	mulps  0x3e0(%rcx),%xmm11
   0x00000000000013ee <+5102>:	addps  %xmm9,%xmm13
   0x0000000000001406 <+5126>:	addps  %xmm8,%xmm15
   0x000000000000144c <+5196>:	movaps 0x3e0(%r8),%xmm10
   0x0000000000001454 <+5204>:	mulps  0x3e0(%rcx),%xmm10
   0x000000000000145c <+5212>:	addps  %xmm7,%xmm11
   0x0000000000001473 <+5235>:	addps  %xmm0,%xmm10

489
490	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_14]);
491	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
492	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000013f2 <+5106>:	movaps 0x330(%r8),%xmm9
   0x00000000000013fa <+5114>:	mulps  0x3f0(%rcx),%xmm9
   0x000000000000140a <+5130>:	movaps 0x370(%r8),%xmm8
   0x0000000000001412 <+5138>:	mulps  0x3f0(%rcx),%xmm8
   0x000000000000141a <+5146>:	addps  %xmm13,%xmm9
   0x0000000000001431 <+5169>:	addps  %xmm15,%xmm8
   0x0000000000001460 <+5216>:	movaps 0x3b0(%r8),%xmm7
   0x0000000000001468 <+5224>:	mulps  0x3f0(%rcx),%xmm7
   0x0000000000001477 <+5239>:	addps  %xmm11,%xmm7
   0x000000000000147b <+5243>:	movaps 0x3f0(%r8),%xmm11
   0x0000000000001483 <+5251>:	mulps  0x3f0(%rcx),%xmm11
   0x000000000000148b <+5259>:	addps  %xmm10,%xmm11

493	      }
494	    }
495
496	    /* Store C(1,4) block. */
497	    for (i = 0; i < 4; i++)
498	    {
499	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x000000000000148f <+5263>:	mulps  %xmm6,%xmm9
   0x00000000000014b5 <+5301>:	mulps  %xmm6,%xmm8
   0x00000000000014cd <+5325>:	mulps  %xmm6,%xmm7
   0x00000000000014e3 <+5347>:	mulps  %xmm6,%xmm11

500	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_14]), C_row[i]);
   0x0000000000001493 <+5267>:	addps  0xc0(%r9),%xmm9
   0x00000000000014b9 <+5305>:	addps  0xd0(%r9),%xmm8
   0x00000000000014d0 <+5328>:	addps  0xe0(%r9),%xmm7
   0x00000000000014e7 <+5351>:	addps  0xf0(%r9),%xmm11

501	      _mm_store_ps(&C[i*4+C_OFFSET_14], C_row[i]);
   0x00000000000014a5 <+5285>:	movaps %xmm9,0xc0(%r9)
   0x00000000000014c1 <+5313>:	movaps %xmm8,0xd0(%r9)
   0x00000000000014d8 <+5336>:	movaps %xmm7,0xe0(%r9)
   0x00000000000014ef <+5359>:	movaps %xmm11,0xf0(%r9)

502	    }
503
504	    /* Reset C(2,1) matrix accumulators */
505	    C_row[0] = _mm_setzero_ps();
   0x00000000000014a1 <+5281>:	xorps  %xmm10,%xmm10

506	    C_row[1] = _mm_setzero_ps();
   0x00000000000014c9 <+5321>:	xorps  %xmm8,%xmm8

507	    C_row[2] = _mm_setzero_ps();
   0x00000000000014e0 <+5344>:	xorps  %xmm7,%xmm7

508	    C_row[3] = _mm_setzero_ps();
   0x00000000000014ad <+5293>:	xorps  %xmm12,%xmm12

509
510	    if (norm[4]*norm[16] >= tolerance &&
   0x000000000000149b <+5275>:	movss  0x28(%rdx,%rsi,1),%xmm0
   0x00000000000014b1 <+5297>:	movaps %xmm5,%xmm9
   0x00000000000014f7 <+5367>:	mulss  %xmm0,%xmm9
   0x00000000000014fc <+5372>:	comiss %xmm1,%xmm9
   0x0000000000001500 <+5376>:	jb     0x1a06 <stream_kernel+6662>

511	        norm[5]*norm[20] >= tolerance &&
   0x0000000000001506 <+5382>:	movss  0x2c(%rdx,%rsi,1),%xmm9
   0x000000000000150d <+5389>:	mulss  0x68(%rdx,%rsi,1),%xmm9
   0x0000000000001514 <+5396>:	comiss %xmm1,%xmm9
   0x0000000000001518 <+5400>:	jb     0x1a06 <stream_kernel+6662>

512	        norm[6]*norm[24] >= tolerance &&
   0x000000000000151e <+5406>:	movss  0x30(%rdx,%rsi,1),%xmm9
   0x0000000000001525 <+5413>:	mulss  0x78(%rdx,%rsi,1),%xmm9
   0x000000000000152c <+5420>:	comiss %xmm1,%xmm9
   0x0000000000001530 <+5424>:	jb     0x1a06 <stream_kernel+6662>

513	        norm[7]*norm[28] >= tolerance)
   0x0000000000001536 <+5430>:	movss  0x34(%rdx,%rsi,1),%xmm9
   0x000000000000153d <+5437>:	mulss  0x88(%rdx,%rsi,1),%xmm9
   0x0000000000001547 <+5447>:	comiss %xmm1,%xmm9
   0x000000000000154b <+5451>:	jb     0x1a06 <stream_kernel+6662>

514	    {
515	      /* A(2,1)*B(1,1) = C(2,1). */
516	      for (i = 0; i < 4; i++)
517	      {
518	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
519	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
520	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001551 <+5457>:	movaps 0x400(%r8),%xmm10
   0x0000000000001571 <+5489>:	movaps 0x440(%r8),%xmm13
   0x0000000000001589 <+5513>:	movaps 0x480(%r8),%xmm8
   0x0000000000001591 <+5521>:	mulps  (%rcx),%xmm10
   0x00000000000015a4 <+5540>:	mulps  (%rcx),%xmm13
   0x00000000000015b1 <+5553>:	mulps  (%rcx),%xmm8
   0x00000000000015f9 <+5625>:	movaps 0x4c0(%r8),%xmm13
   0x0000000000001601 <+5633>:	mulps  (%rcx),%xmm13

521
522	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
523	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
524	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001559 <+5465>:	movaps 0x410(%r8),%xmm11
   0x0000000000001579 <+5497>:	movaps 0x450(%r8),%xmm14
   0x0000000000001595 <+5525>:	mulps  0x10(%rcx),%xmm11
   0x00000000000015a8 <+5544>:	mulps  0x10(%rcx),%xmm14
   0x00000000000015b5 <+5557>:	movaps 0x4d0(%r8),%xmm15
   0x00000000000015bd <+5565>:	addps  %xmm10,%xmm11
   0x00000000000015c1 <+5569>:	mulps  0x10(%rcx),%xmm15
   0x00000000000015d2 <+5586>:	movaps 0x490(%r8),%xmm11
   0x00000000000015df <+5599>:	mulps  0x10(%rcx),%xmm11
   0x00000000000015f5 <+5621>:	addps  %xmm13,%xmm14
   0x000000000000162a <+5674>:	addps  %xmm8,%xmm11
   0x0000000000001663 <+5731>:	addps  %xmm13,%xmm15

525
526	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
527	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
528	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001561 <+5473>:	movaps 0x420(%r8),%xmm12
   0x000000000000159a <+5530>:	mulps  0x20(%rcx),%xmm12
   0x00000000000015c6 <+5574>:	movaps 0x460(%r8),%xmm10
   0x00000000000015ce <+5582>:	addps  %xmm11,%xmm12
   0x00000000000015da <+5594>:	mulps  0x20(%rcx),%xmm10
   0x00000000000015e8 <+5608>:	movaps 0x4a0(%r8),%xmm12
   0x00000000000015f0 <+5616>:	mulps  0x20(%rcx),%xmm12
   0x0000000000001605 <+5637>:	addps  %xmm14,%xmm10
   0x000000000000161d <+5661>:	movaps 0x4e0(%r8),%xmm10
   0x0000000000001625 <+5669>:	mulps  0x20(%rcx),%xmm10
   0x000000000000163b <+5691>:	addps  %xmm11,%xmm12
   0x0000000000001677 <+5751>:	addps  %xmm15,%xmm10

529
530	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
531	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
532	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001569 <+5481>:	movaps 0x430(%r8),%xmm9
   0x0000000000001581 <+5505>:	movaps 0x470(%r8),%xmm7
   0x000000000000159f <+5535>:	mulps  0x30(%rcx),%xmm9
   0x00000000000015ad <+5549>:	mulps  0x30(%rcx),%xmm7
   0x00000000000015e4 <+5604>:	addps  %xmm12,%xmm9
   0x0000000000001619 <+5657>:	addps  %xmm10,%xmm7
   0x000000000000162e <+5678>:	movaps 0x4b0(%r8),%xmm8
   0x0000000000001636 <+5686>:	mulps  0x30(%rcx),%xmm8
   0x000000000000164f <+5711>:	addps  %xmm12,%xmm8
   0x000000000000167b <+5755>:	movaps 0x4f0(%r8),%xmm15
   0x0000000000001683 <+5763>:	mulps  0x30(%rcx),%xmm15
   0x000000000000169c <+5788>:	addps  %xmm10,%xmm15

533	      }
534
535	      /* A(2,2)*B(2,1) = C(2,1). */
536	      for (i = 0; i < 4; i++)
537	      {
538	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
539	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
540	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000163f <+5695>:	movaps 0x500(%r8),%xmm11
   0x0000000000001647 <+5703>:	mulps  0x100(%rcx),%xmm11
   0x0000000000001688 <+5768>:	addps  %xmm9,%xmm11
   0x00000000000016b4 <+5812>:	movaps 0x540(%r8),%xmm11
   0x00000000000016bc <+5820>:	mulps  0x100(%rcx),%xmm11
   0x00000000000016c8 <+5832>:	movaps 0x580(%r8),%xmm9
   0x00000000000016d0 <+5840>:	mulps  0x100(%rcx),%xmm9
   0x00000000000016ec <+5868>:	addps  %xmm7,%xmm11
   0x00000000000016ff <+5887>:	addps  %xmm8,%xmm9
   0x000000000000173e <+5950>:	movaps 0x5c0(%r8),%xmm10
   0x0000000000001746 <+5958>:	mulps  0x100(%rcx),%xmm10
   0x000000000000178a <+6026>:	addps  %xmm15,%xmm10

541
542	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
543	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
544	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000168c <+5772>:	movaps 0x510(%r8),%xmm9
   0x0000000000001694 <+5780>:	mulps  0x110(%rcx),%xmm9
   0x00000000000016b0 <+5808>:	addps  %xmm11,%xmm9
   0x00000000000016f0 <+5872>:	movaps 0x550(%r8),%xmm7
   0x00000000000016f8 <+5880>:	mulps  0x110(%rcx),%xmm7
   0x0000000000001703 <+5891>:	movaps 0x590(%r8),%xmm8
   0x000000000000170b <+5899>:	mulps  0x110(%rcx),%xmm8
   0x0000000000001713 <+5907>:	addps  %xmm11,%xmm7
   0x000000000000174e <+5966>:	addps  %xmm9,%xmm8
   0x000000000000178e <+6030>:	movaps 0x5d0(%r8),%xmm15
   0x0000000000001796 <+6038>:	mulps  0x110(%rcx),%xmm15
   0x00000000000017b2 <+6066>:	addps  %xmm10,%xmm15

545
546	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
547	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
548	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000016a0 <+5792>:	movaps 0x520(%r8),%xmm10
   0x00000000000016a8 <+5800>:	mulps  0x120(%rcx),%xmm10
   0x00000000000016c4 <+5828>:	addps  %xmm9,%xmm10
   0x00000000000016dc <+5852>:	movaps 0x560(%r8),%xmm10
   0x00000000000016e4 <+5860>:	mulps  0x120(%rcx),%xmm10
   0x0000000000001717 <+5911>:	movaps 0x5a0(%r8),%xmm11
   0x000000000000171f <+5919>:	mulps  0x120(%rcx),%xmm11
   0x0000000000001727 <+5927>:	addps  %xmm7,%xmm10
   0x0000000000001762 <+5986>:	addps  %xmm8,%xmm11
   0x00000000000017b6 <+6070>:	movaps 0x5e0(%r8),%xmm10
   0x00000000000017be <+6078>:	mulps  0x120(%rcx),%xmm10
   0x00000000000017da <+6106>:	addps  %xmm15,%xmm10

549
550	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
551	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
552	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001609 <+5641>:	movaps 0x530(%r8),%xmm14
   0x0000000000001611 <+5649>:	mulps  0x130(%rcx),%xmm14
   0x0000000000001653 <+5715>:	movaps 0x5b0(%r8),%xmm12
   0x000000000000165b <+5723>:	mulps  0x130(%rcx),%xmm12
   0x0000000000001667 <+5735>:	movaps 0x570(%r8),%xmm13
   0x000000000000166f <+5743>:	mulps  0x130(%rcx),%xmm13
   0x00000000000016d8 <+5848>:	addps  %xmm10,%xmm14
   0x000000000000173a <+5946>:	addps  %xmm10,%xmm13
   0x0000000000001776 <+6006>:	addps  %xmm11,%xmm12
   0x000000000000177a <+6010>:	movaps 0x5f0(%r8),%xmm11
   0x0000000000001782 <+6018>:	mulps  0x130(%rcx),%xmm11
   0x00000000000017ee <+6126>:	addps  %xmm10,%xmm11

553	      }
554
555	      /* A(2,3)*B(3,1) = C(2,1). */
556	      for (i = 0; i < 4; i++)
557	      {
558	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
559	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
560	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000172b <+5931>:	movaps 0x640(%r8),%xmm7
   0x0000000000001733 <+5939>:	mulps  0x200(%rcx),%xmm7
   0x0000000000001752 <+5970>:	movaps 0x680(%r8),%xmm9
   0x000000000000175a <+5978>:	mulps  0x200(%rcx),%xmm9
   0x0000000000001766 <+5990>:	movaps 0x600(%r8),%xmm8
   0x000000000000176e <+5998>:	mulps  0x200(%rcx),%xmm8
   0x000000000000179e <+6046>:	addps  %xmm14,%xmm8
   0x000000000000182a <+6186>:	addps  %xmm13,%xmm7
   0x000000000000183e <+6206>:	addps  %xmm12,%xmm9
   0x000000000000187d <+6269>:	movaps 0x6c0(%r8),%xmm13
   0x0000000000001885 <+6277>:	mulps  0x200(%rcx),%xmm13
   0x00000000000018c8 <+6344>:	addps  %xmm11,%xmm13

561
562	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
563	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
564	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000017a2 <+6050>:	movaps 0x610(%r8),%xmm14
   0x00000000000017aa <+6058>:	mulps  0x210(%rcx),%xmm14
   0x00000000000017c6 <+6086>:	addps  %xmm8,%xmm14
   0x000000000000182e <+6190>:	movaps 0x650(%r8),%xmm13
   0x0000000000001836 <+6198>:	mulps  0x210(%rcx),%xmm13
   0x0000000000001842 <+6210>:	movaps 0x690(%r8),%xmm12
   0x000000000000184a <+6218>:	mulps  0x210(%rcx),%xmm12
   0x0000000000001852 <+6226>:	addps  %xmm7,%xmm13
   0x0000000000001865 <+6245>:	addps  %xmm9,%xmm12
   0x00000000000018cc <+6348>:	movaps 0x6d0(%r8),%xmm11
   0x00000000000018d4 <+6356>:	mulps  0x210(%rcx),%xmm11
   0x00000000000018f0 <+6384>:	addps  %xmm13,%xmm11

565
566	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
567	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
568	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000017ca <+6090>:	movaps 0x620(%r8),%xmm8
   0x00000000000017d2 <+6098>:	mulps  0x220(%rcx),%xmm8
   0x0000000000001802 <+6146>:	addps  %xmm14,%xmm8
   0x0000000000001806 <+6150>:	movaps 0x6e0(%r8),%xmm14
   0x000000000000180e <+6158>:	mulps  0x220(%rcx),%xmm14
   0x0000000000001856 <+6230>:	movaps 0x660(%r8),%xmm7
   0x000000000000185e <+6238>:	mulps  0x220(%rcx),%xmm7
   0x0000000000001869 <+6249>:	movaps 0x6a0(%r8),%xmm9
   0x0000000000001871 <+6257>:	mulps  0x220(%rcx),%xmm9
   0x0000000000001879 <+6265>:	addps  %xmm13,%xmm7
   0x00000000000018a0 <+6304>:	addps  %xmm12,%xmm9
   0x0000000000001904 <+6404>:	addps  %xmm11,%xmm14

569
570	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
571	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
572	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000017f2 <+6130>:	movaps 0x630(%r8),%xmm10
   0x00000000000017fa <+6138>:	mulps  0x230(%rcx),%xmm10
   0x0000000000001816 <+6166>:	addps  %xmm8,%xmm10
   0x000000000000181a <+6170>:	movaps 0x670(%r8),%xmm8
   0x0000000000001822 <+6178>:	mulps  0x230(%rcx),%xmm8
   0x000000000000188d <+6285>:	addps  %xmm7,%xmm8
   0x0000000000001891 <+6289>:	movaps 0x6b0(%r8),%xmm7
   0x0000000000001899 <+6297>:	mulps  0x230(%rcx),%xmm7
   0x00000000000018b4 <+6324>:	addps  %xmm9,%xmm7
   0x00000000000018b8 <+6328>:	movaps 0x6f0(%r8),%xmm9
   0x00000000000018c0 <+6336>:	mulps  0x230(%rcx),%xmm9
   0x0000000000001918 <+6424>:	addps  %xmm14,%xmm9

573	      }
574
575	      /* A(2,4)*B(4,1) = C(2,1). */
576	      for (i = 0; i < 4; i++)
577	      {
578	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
579	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
580	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000017de <+6110>:	movaps 0x700(%r8),%xmm15
   0x00000000000017e6 <+6118>:	mulps  0x300(%rcx),%xmm15
   0x00000000000018a4 <+6308>:	movaps 0x740(%r8),%xmm12
   0x00000000000018ac <+6316>:	mulps  0x300(%rcx),%xmm12
   0x00000000000018dc <+6364>:	addps  %xmm10,%xmm15
   0x000000000000191c <+6428>:	movaps 0x780(%r8),%xmm14
   0x0000000000001924 <+6436>:	mulps  0x300(%rcx),%xmm14
   0x0000000000001944 <+6468>:	addps  %xmm8,%xmm12
   0x0000000000001970 <+6512>:	movaps 0x7c0(%r8),%xmm12
   0x0000000000001978 <+6520>:	mulps  0x300(%rcx),%xmm12
   0x0000000000001994 <+6548>:	addps  %xmm7,%xmm14
   0x00000000000019ab <+6571>:	addps  %xmm9,%xmm12

581
582	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
583	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
584	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000018e0 <+6368>:	movaps 0x710(%r8),%xmm10
   0x00000000000018e8 <+6376>:	mulps  0x310(%rcx),%xmm10
   0x000000000000192c <+6444>:	addps  %xmm15,%xmm10
   0x0000000000001948 <+6472>:	movaps 0x750(%r8),%xmm8
   0x0000000000001950 <+6480>:	mulps  0x310(%rcx),%xmm8
   0x000000000000196c <+6508>:	addps  %xmm12,%xmm8
   0x0000000000001998 <+6552>:	movaps 0x790(%r8),%xmm7
   0x00000000000019a0 <+6560>:	mulps  0x310(%rcx),%xmm7
   0x00000000000019af <+6575>:	movaps 0x7d0(%r8),%xmm9
   0x00000000000019b7 <+6583>:	mulps  0x310(%rcx),%xmm9
   0x00000000000019bf <+6591>:	addps  %xmm14,%xmm7
   0x00000000000019d6 <+6614>:	addps  %xmm12,%xmm9

585
586	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
587	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
588	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000018f4 <+6388>:	movaps 0x760(%r8),%xmm13
   0x00000000000018fc <+6396>:	mulps  0x320(%rcx),%xmm13
   0x0000000000001908 <+6408>:	movaps 0x720(%r8),%xmm11
   0x0000000000001910 <+6416>:	mulps  0x320(%rcx),%xmm11
   0x0000000000001930 <+6448>:	addps  %xmm10,%xmm11
   0x000000000000195c <+6492>:	movaps 0x7a0(%r8),%xmm11
   0x0000000000001964 <+6500>:	mulps  0x320(%rcx),%xmm11
   0x0000000000001980 <+6528>:	addps  %xmm8,%xmm13
   0x00000000000019c3 <+6595>:	addps  %xmm7,%xmm11
   0x00000000000019ee <+6638>:	movaps 0x7e0(%r8),%xmm11
   0x00000000000019f6 <+6646>:	mulps  0x320(%rcx),%xmm11
   0x00000000000019fe <+6654>:	addps  %xmm9,%xmm11

589
590	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
591	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
592	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001934 <+6452>:	movaps 0x730(%r8),%xmm10
   0x000000000000193c <+6460>:	mulps  0x330(%rcx),%xmm10
   0x0000000000001958 <+6488>:	addps  %xmm11,%xmm10
   0x0000000000001984 <+6532>:	movaps 0x770(%r8),%xmm8
   0x000000000000198c <+6540>:	mulps  0x330(%rcx),%xmm8
   0x00000000000019a7 <+6567>:	addps  %xmm13,%xmm8
   0x00000000000019c7 <+6599>:	movaps 0x7b0(%r8),%xmm7
   0x00000000000019cf <+6607>:	mulps  0x330(%rcx),%xmm7
   0x00000000000019da <+6618>:	movaps 0x7f0(%r8),%xmm12
   0x00000000000019e2 <+6626>:	mulps  0x330(%rcx),%xmm12
   0x00000000000019ea <+6634>:	addps  %xmm11,%xmm7
   0x0000000000001a02 <+6658>:	addps  %xmm11,%xmm12

593	      }
594	    }
595
596	    /* Store C(2,1) block. */
597	    for (i = 0; i < 4; i++)
598	    {
599	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000001a06 <+6662>:	mulps  %xmm6,%xmm10
   0x0000000000001a22 <+6690>:	mulps  %xmm6,%xmm8
   0x0000000000001a3a <+6714>:	mulps  %xmm6,%xmm7
   0x0000000000001a50 <+6736>:	mulps  %xmm6,%xmm12

600	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_21]), C_row[i]);
   0x0000000000001a0a <+6666>:	addps  0x100(%r9),%xmm10
   0x0000000000001a26 <+6694>:	addps  0x110(%r9),%xmm8
   0x0000000000001a3d <+6717>:	addps  0x120(%r9),%xmm7
   0x0000000000001a54 <+6740>:	addps  0x130(%r9),%xmm12

601	      _mm_store_ps(&C[i*4+C_OFFSET_21], C_row[i]);
   0x0000000000001a16 <+6678>:	movaps %xmm10,0x100(%r9)
   0x0000000000001a2e <+6702>:	movaps %xmm8,0x110(%r9)
   0x0000000000001a45 <+6725>:	movaps %xmm7,0x120(%r9)
   0x0000000000001a5c <+6748>:	movaps %xmm12,0x130(%r9)

602	    }
603
604	    /* Reset C(2,2) matrix accumulators */
605	    C_row[0] = _mm_setzero_ps();
   0x0000000000001a1e <+6686>:	xorps  %xmm10,%xmm10

606	    C_row[1] = _mm_setzero_ps();
   0x0000000000001a36 <+6710>:	xorps  %xmm8,%xmm8

607	    C_row[2] = _mm_setzero_ps();
   0x0000000000001a4d <+6733>:	xorps  %xmm7,%xmm7

608	    C_row[3] = _mm_setzero_ps();
   0x0000000000001a64 <+6756>:	xorps  %xmm12,%xmm12

609
610	    if (norm[4]*norm[17] >= tolerance &&
   0x0000000000001a12 <+6674>:	movaps %xmm4,%xmm9
   0x0000000000001a68 <+6760>:	mulss  %xmm0,%xmm9
   0x0000000000001a6d <+6765>:	comiss %xmm1,%xmm9
   0x0000000000001a71 <+6769>:	jb     0x1f7b <stream_kernel+8059>

611	        norm[5]*norm[21] >= tolerance &&
   0x0000000000001a77 <+6775>:	movss  0x2c(%rdx,%rsi,1),%xmm9
   0x0000000000001a7e <+6782>:	mulss  0x6c(%rdx,%rsi,1),%xmm9
   0x0000000000001a85 <+6789>:	comiss %xmm1,%xmm9
   0x0000000000001a89 <+6793>:	jb     0x1f7b <stream_kernel+8059>

612	        norm[6]*norm[25] >= tolerance &&
   0x0000000000001a8f <+6799>:	movss  0x30(%rdx,%rsi,1),%xmm9
   0x0000000000001a96 <+6806>:	mulss  0x7c(%rdx,%rsi,1),%xmm9
   0x0000000000001a9d <+6813>:	comiss %xmm1,%xmm9
   0x0000000000001aa1 <+6817>:	jb     0x1f7b <stream_kernel+8059>

613	        norm[7]*norm[29] >= tolerance)
   0x0000000000001aa7 <+6823>:	movss  0x34(%rdx,%rsi,1),%xmm9
   0x0000000000001aae <+6830>:	mulss  0x8c(%rdx,%rsi,1),%xmm9
   0x0000000000001ab8 <+6840>:	comiss %xmm1,%xmm9
   0x0000000000001abc <+6844>:	jb     0x1f7b <stream_kernel+8059>

614	    {
615	      /* A(2,1)*B(1,2) = C(2,2). */
616	      for (i = 0; i < 4; i++)
617	      {
618	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
619	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
620	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001ac2 <+6850>:	movaps 0x400(%r8),%xmm10
   0x0000000000001ae2 <+6882>:	movaps 0x440(%r8),%xmm13
   0x0000000000001afa <+6906>:	movaps 0x480(%r8),%xmm8
   0x0000000000001b02 <+6914>:	mulps  0x40(%rcx),%xmm10
   0x0000000000001b16 <+6934>:	mulps  0x40(%rcx),%xmm13
   0x0000000000001b24 <+6948>:	mulps  0x40(%rcx),%xmm8
   0x0000000000001b6d <+7021>:	movaps 0x4c0(%r8),%xmm13
   0x0000000000001b75 <+7029>:	mulps  0x40(%rcx),%xmm13

621
622	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
623	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
624	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001aca <+6858>:	movaps 0x410(%r8),%xmm11
   0x0000000000001aea <+6890>:	movaps 0x450(%r8),%xmm14
   0x0000000000001b07 <+6919>:	mulps  0x50(%rcx),%xmm11
   0x0000000000001b1b <+6939>:	mulps  0x50(%rcx),%xmm14
   0x0000000000001b29 <+6953>:	movaps 0x4d0(%r8),%xmm15
   0x0000000000001b31 <+6961>:	addps  %xmm10,%xmm11
   0x0000000000001b35 <+6965>:	mulps  0x50(%rcx),%xmm15
   0x0000000000001b46 <+6982>:	movaps 0x490(%r8),%xmm11
   0x0000000000001b53 <+6995>:	mulps  0x50(%rcx),%xmm11
   0x0000000000001b69 <+7017>:	addps  %xmm13,%xmm14
   0x0000000000001b9f <+7071>:	addps  %xmm8,%xmm11
   0x0000000000001bd8 <+7128>:	addps  %xmm13,%xmm15

625
626	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
627	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
628	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001ad2 <+6866>:	movaps 0x420(%r8),%xmm12
   0x0000000000001b0c <+6924>:	mulps  0x60(%rcx),%xmm12
   0x0000000000001b3a <+6970>:	movaps 0x460(%r8),%xmm10
   0x0000000000001b42 <+6978>:	addps  %xmm11,%xmm12
   0x0000000000001b4e <+6990>:	mulps  0x60(%rcx),%xmm10
   0x0000000000001b5c <+7004>:	movaps 0x4a0(%r8),%xmm12
   0x0000000000001b64 <+7012>:	mulps  0x60(%rcx),%xmm12
   0x0000000000001b7a <+7034>:	addps  %xmm14,%xmm10
   0x0000000000001b92 <+7058>:	movaps 0x4e0(%r8),%xmm10
   0x0000000000001b9a <+7066>:	mulps  0x60(%rcx),%xmm10
   0x0000000000001bb0 <+7088>:	addps  %xmm11,%xmm12
   0x0000000000001bec <+7148>:	addps  %xmm15,%xmm10

629
630	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
631	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
632	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001ada <+6874>:	movaps 0x430(%r8),%xmm9
   0x0000000000001af2 <+6898>:	movaps 0x470(%r8),%xmm7
   0x0000000000001b11 <+6929>:	mulps  0x70(%rcx),%xmm9
   0x0000000000001b20 <+6944>:	mulps  0x70(%rcx),%xmm7
   0x0000000000001b58 <+7000>:	addps  %xmm12,%xmm9
   0x0000000000001b8e <+7054>:	addps  %xmm10,%xmm7
   0x0000000000001ba3 <+7075>:	movaps 0x4b0(%r8),%xmm8
   0x0000000000001bab <+7083>:	mulps  0x70(%rcx),%xmm8
   0x0000000000001bc4 <+7108>:	addps  %xmm12,%xmm8
   0x0000000000001bf0 <+7152>:	movaps 0x4f0(%r8),%xmm15
   0x0000000000001bf8 <+7160>:	mulps  0x70(%rcx),%xmm15
   0x0000000000001c11 <+7185>:	addps  %xmm10,%xmm15

633	      }
634
635	      /* A(2,2)*B(2,2) = C(2,2). */
636	      for (i = 0; i < 4; i++)
637	      {
638	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
639	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
640	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001bb4 <+7092>:	movaps 0x500(%r8),%xmm11
   0x0000000000001bbc <+7100>:	mulps  0x140(%rcx),%xmm11
   0x0000000000001bfd <+7165>:	addps  %xmm9,%xmm11
   0x0000000000001c29 <+7209>:	movaps 0x540(%r8),%xmm11
   0x0000000000001c31 <+7217>:	mulps  0x140(%rcx),%xmm11
   0x0000000000001c3d <+7229>:	movaps 0x580(%r8),%xmm9
   0x0000000000001c45 <+7237>:	mulps  0x140(%rcx),%xmm9
   0x0000000000001c61 <+7265>:	addps  %xmm7,%xmm11
   0x0000000000001c74 <+7284>:	addps  %xmm8,%xmm9
   0x0000000000001cb3 <+7347>:	movaps 0x5c0(%r8),%xmm10
   0x0000000000001cbb <+7355>:	mulps  0x140(%rcx),%xmm10
   0x0000000000001cff <+7423>:	addps  %xmm15,%xmm10

641
642	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
643	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
644	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c01 <+7169>:	movaps 0x510(%r8),%xmm9
   0x0000000000001c09 <+7177>:	mulps  0x150(%rcx),%xmm9
   0x0000000000001c25 <+7205>:	addps  %xmm11,%xmm9
   0x0000000000001c65 <+7269>:	movaps 0x550(%r8),%xmm7
   0x0000000000001c6d <+7277>:	mulps  0x150(%rcx),%xmm7
   0x0000000000001c78 <+7288>:	movaps 0x590(%r8),%xmm8
   0x0000000000001c80 <+7296>:	mulps  0x150(%rcx),%xmm8
   0x0000000000001c88 <+7304>:	addps  %xmm11,%xmm7
   0x0000000000001cc3 <+7363>:	addps  %xmm9,%xmm8
   0x0000000000001d03 <+7427>:	movaps 0x5d0(%r8),%xmm15
   0x0000000000001d0b <+7435>:	mulps  0x150(%rcx),%xmm15
   0x0000000000001d27 <+7463>:	addps  %xmm10,%xmm15

645
646	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
647	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
648	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001c15 <+7189>:	movaps 0x520(%r8),%xmm10
   0x0000000000001c1d <+7197>:	mulps  0x160(%rcx),%xmm10
   0x0000000000001c39 <+7225>:	addps  %xmm9,%xmm10
   0x0000000000001c51 <+7249>:	movaps 0x560(%r8),%xmm10
   0x0000000000001c59 <+7257>:	mulps  0x160(%rcx),%xmm10
   0x0000000000001c8c <+7308>:	movaps 0x5a0(%r8),%xmm11
   0x0000000000001c94 <+7316>:	mulps  0x160(%rcx),%xmm11
   0x0000000000001c9c <+7324>:	addps  %xmm7,%xmm10
   0x0000000000001cd7 <+7383>:	addps  %xmm8,%xmm11
   0x0000000000001d2b <+7467>:	movaps 0x5e0(%r8),%xmm10
   0x0000000000001d33 <+7475>:	mulps  0x160(%rcx),%xmm10
   0x0000000000001d4f <+7503>:	addps  %xmm15,%xmm10

649
650	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
651	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
652	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001b7e <+7038>:	movaps 0x530(%r8),%xmm14
   0x0000000000001b86 <+7046>:	mulps  0x170(%rcx),%xmm14
   0x0000000000001bc8 <+7112>:	movaps 0x5b0(%r8),%xmm12
   0x0000000000001bd0 <+7120>:	mulps  0x170(%rcx),%xmm12
   0x0000000000001bdc <+7132>:	movaps 0x570(%r8),%xmm13
   0x0000000000001be4 <+7140>:	mulps  0x170(%rcx),%xmm13
   0x0000000000001c4d <+7245>:	addps  %xmm10,%xmm14
   0x0000000000001caf <+7343>:	addps  %xmm10,%xmm13
   0x0000000000001ceb <+7403>:	addps  %xmm11,%xmm12
   0x0000000000001cef <+7407>:	movaps 0x5f0(%r8),%xmm11
   0x0000000000001cf7 <+7415>:	mulps  0x170(%rcx),%xmm11
   0x0000000000001d63 <+7523>:	addps  %xmm10,%xmm11

653	      }
654
655	      /* A(2,3)*B(3,2) = C(2,2). */
656	      for (i = 0; i < 4; i++)
657	      {
658	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
659	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
660	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001ca0 <+7328>:	movaps 0x640(%r8),%xmm7
   0x0000000000001ca8 <+7336>:	mulps  0x240(%rcx),%xmm7
   0x0000000000001cc7 <+7367>:	movaps 0x680(%r8),%xmm9
   0x0000000000001ccf <+7375>:	mulps  0x240(%rcx),%xmm9
   0x0000000000001cdb <+7387>:	movaps 0x600(%r8),%xmm8
   0x0000000000001ce3 <+7395>:	mulps  0x240(%rcx),%xmm8
   0x0000000000001d13 <+7443>:	addps  %xmm14,%xmm8
   0x0000000000001d9f <+7583>:	addps  %xmm13,%xmm7
   0x0000000000001db3 <+7603>:	addps  %xmm12,%xmm9
   0x0000000000001df2 <+7666>:	movaps 0x6c0(%r8),%xmm13
   0x0000000000001dfa <+7674>:	mulps  0x240(%rcx),%xmm13
   0x0000000000001e3d <+7741>:	addps  %xmm11,%xmm13

661
662	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
663	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
664	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001d17 <+7447>:	movaps 0x610(%r8),%xmm14
   0x0000000000001d1f <+7455>:	mulps  0x250(%rcx),%xmm14
   0x0000000000001d3b <+7483>:	addps  %xmm8,%xmm14
   0x0000000000001da3 <+7587>:	movaps 0x650(%r8),%xmm13
   0x0000000000001dab <+7595>:	mulps  0x250(%rcx),%xmm13
   0x0000000000001db7 <+7607>:	movaps 0x690(%r8),%xmm12
   0x0000000000001dbf <+7615>:	mulps  0x250(%rcx),%xmm12
   0x0000000000001dc7 <+7623>:	addps  %xmm7,%xmm13
   0x0000000000001dda <+7642>:	addps  %xmm9,%xmm12
   0x0000000000001e41 <+7745>:	movaps 0x6d0(%r8),%xmm11
   0x0000000000001e49 <+7753>:	mulps  0x250(%rcx),%xmm11
   0x0000000000001e65 <+7781>:	addps  %xmm13,%xmm11

665
666	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
667	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
668	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001d3f <+7487>:	movaps 0x620(%r8),%xmm8
   0x0000000000001d47 <+7495>:	mulps  0x260(%rcx),%xmm8
   0x0000000000001d77 <+7543>:	addps  %xmm14,%xmm8
   0x0000000000001d7b <+7547>:	movaps 0x6e0(%r8),%xmm14
   0x0000000000001d83 <+7555>:	mulps  0x260(%rcx),%xmm14
   0x0000000000001dcb <+7627>:	movaps 0x660(%r8),%xmm7
   0x0000000000001dd3 <+7635>:	mulps  0x260(%rcx),%xmm7
   0x0000000000001dde <+7646>:	movaps 0x6a0(%r8),%xmm9
   0x0000000000001de6 <+7654>:	mulps  0x260(%rcx),%xmm9
   0x0000000000001dee <+7662>:	addps  %xmm13,%xmm7
   0x0000000000001e15 <+7701>:	addps  %xmm12,%xmm9
   0x0000000000001e79 <+7801>:	addps  %xmm11,%xmm14

669
670	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
671	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
672	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001d67 <+7527>:	movaps 0x630(%r8),%xmm10
   0x0000000000001d6f <+7535>:	mulps  0x270(%rcx),%xmm10
   0x0000000000001d8b <+7563>:	addps  %xmm8,%xmm10
   0x0000000000001d8f <+7567>:	movaps 0x670(%r8),%xmm8
   0x0000000000001d97 <+7575>:	mulps  0x270(%rcx),%xmm8
   0x0000000000001e02 <+7682>:	addps  %xmm7,%xmm8
   0x0000000000001e06 <+7686>:	movaps 0x6b0(%r8),%xmm7
   0x0000000000001e0e <+7694>:	mulps  0x270(%rcx),%xmm7
   0x0000000000001e29 <+7721>:	addps  %xmm9,%xmm7
   0x0000000000001e2d <+7725>:	movaps 0x6f0(%r8),%xmm9
   0x0000000000001e35 <+7733>:	mulps  0x270(%rcx),%xmm9
   0x0000000000001e8d <+7821>:	addps  %xmm14,%xmm9

673	      }
674
675	      /* A(2,4)*B(4,2) = C(2,2). */
676	      for (i = 0; i < 4; i++)
677	      {
678	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
679	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
680	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001d53 <+7507>:	movaps 0x700(%r8),%xmm15
   0x0000000000001d5b <+7515>:	mulps  0x340(%rcx),%xmm15
   0x0000000000001e19 <+7705>:	movaps 0x740(%r8),%xmm12
   0x0000000000001e21 <+7713>:	mulps  0x340(%rcx),%xmm12
   0x0000000000001e51 <+7761>:	addps  %xmm10,%xmm15
   0x0000000000001e91 <+7825>:	movaps 0x780(%r8),%xmm14
   0x0000000000001e99 <+7833>:	mulps  0x340(%rcx),%xmm14
   0x0000000000001eb9 <+7865>:	addps  %xmm8,%xmm12
   0x0000000000001ee5 <+7909>:	movaps 0x7c0(%r8),%xmm12
   0x0000000000001eed <+7917>:	mulps  0x340(%rcx),%xmm12
   0x0000000000001f09 <+7945>:	addps  %xmm7,%xmm14
   0x0000000000001f20 <+7968>:	addps  %xmm9,%xmm12

681
682	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
683	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
684	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001e55 <+7765>:	movaps 0x710(%r8),%xmm10
   0x0000000000001e5d <+7773>:	mulps  0x350(%rcx),%xmm10
   0x0000000000001ea1 <+7841>:	addps  %xmm15,%xmm10
   0x0000000000001ebd <+7869>:	movaps 0x750(%r8),%xmm8
   0x0000000000001ec5 <+7877>:	mulps  0x350(%rcx),%xmm8
   0x0000000000001ee1 <+7905>:	addps  %xmm12,%xmm8
   0x0000000000001f0d <+7949>:	movaps 0x790(%r8),%xmm7
   0x0000000000001f15 <+7957>:	mulps  0x350(%rcx),%xmm7
   0x0000000000001f24 <+7972>:	movaps 0x7d0(%r8),%xmm9
   0x0000000000001f2c <+7980>:	mulps  0x350(%rcx),%xmm9
   0x0000000000001f34 <+7988>:	addps  %xmm14,%xmm7
   0x0000000000001f4b <+8011>:	addps  %xmm12,%xmm9

685
686	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
687	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
688	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001e69 <+7785>:	movaps 0x760(%r8),%xmm13
   0x0000000000001e71 <+7793>:	mulps  0x360(%rcx),%xmm13
   0x0000000000001e7d <+7805>:	movaps 0x720(%r8),%xmm11
   0x0000000000001e85 <+7813>:	mulps  0x360(%rcx),%xmm11
   0x0000000000001ea5 <+7845>:	addps  %xmm10,%xmm11
   0x0000000000001ed1 <+7889>:	movaps 0x7a0(%r8),%xmm11
   0x0000000000001ed9 <+7897>:	mulps  0x360(%rcx),%xmm11
   0x0000000000001ef5 <+7925>:	addps  %xmm8,%xmm13
   0x0000000000001f38 <+7992>:	addps  %xmm7,%xmm11
   0x0000000000001f63 <+8035>:	movaps 0x7e0(%r8),%xmm11
   0x0000000000001f6b <+8043>:	mulps  0x360(%rcx),%xmm11
   0x0000000000001f73 <+8051>:	addps  %xmm9,%xmm11

689
690	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
691	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
692	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000001ea9 <+7849>:	movaps 0x730(%r8),%xmm10
   0x0000000000001eb1 <+7857>:	mulps  0x370(%rcx),%xmm10
   0x0000000000001ecd <+7885>:	addps  %xmm11,%xmm10
   0x0000000000001ef9 <+7929>:	movaps 0x770(%r8),%xmm8
   0x0000000000001f01 <+7937>:	mulps  0x370(%rcx),%xmm8
   0x0000000000001f1c <+7964>:	addps  %xmm13,%xmm8
   0x0000000000001f3c <+7996>:	movaps 0x7b0(%r8),%xmm7
   0x0000000000001f44 <+8004>:	mulps  0x370(%rcx),%xmm7
   0x0000000000001f4f <+8015>:	movaps 0x7f0(%r8),%xmm12
   0x0000000000001f57 <+8023>:	mulps  0x370(%rcx),%xmm12
   0x0000000000001f5f <+8031>:	addps  %xmm11,%xmm7
   0x0000000000001f77 <+8055>:	addps  %xmm11,%xmm12

693	      }
694	    }
695
696	    /* Store C(2,2) block. */
697	    for (i = 0; i < 4; i++)
698	    {
699	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000001f7b <+8059>:	mulps  %xmm6,%xmm10
   0x0000000000001f97 <+8087>:	mulps  %xmm6,%xmm8
   0x0000000000001faf <+8111>:	mulps  %xmm6,%xmm7
   0x0000000000001fc5 <+8133>:	mulps  %xmm6,%xmm12

700	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_22]), C_row[i]);
   0x0000000000001f7f <+8063>:	addps  0x140(%r9),%xmm10
   0x0000000000001f9b <+8091>:	addps  0x150(%r9),%xmm8
   0x0000000000001fb2 <+8114>:	addps  0x160(%r9),%xmm7
   0x0000000000001fc9 <+8137>:	addps  0x170(%r9),%xmm12

701	      _mm_store_ps(&C[i*4+C_OFFSET_22], C_row[i]);
   0x0000000000001f8b <+8075>:	movaps %xmm10,0x140(%r9)
   0x0000000000001fa3 <+8099>:	movaps %xmm8,0x150(%r9)
   0x0000000000001fba <+8122>:	movaps %xmm7,0x160(%r9)
   0x0000000000001fd1 <+8145>:	movaps %xmm12,0x170(%r9)

702	    }
703
704	    /* Reset C(2,3) matrix accumulators */
705	    C_row[0] = _mm_setzero_ps();
   0x0000000000001f93 <+8083>:	xorps  %xmm10,%xmm10

706	    C_row[1] = _mm_setzero_ps();
   0x0000000000001fab <+8107>:	xorps  %xmm8,%xmm8

707	    C_row[2] = _mm_setzero_ps();
   0x0000000000001fc2 <+8130>:	xorps  %xmm7,%xmm7

708	    C_row[3] = _mm_setzero_ps();
   0x0000000000001fd9 <+8153>:	xorps  %xmm12,%xmm12

709
710	    if (norm[4]*norm[18] >= tolerance &&
   0x0000000000001f87 <+8071>:	movaps %xmm3,%xmm9
   0x0000000000001fdd <+8157>:	mulss  %xmm0,%xmm9
   0x0000000000001fe2 <+8162>:	comiss %xmm1,%xmm9
   0x0000000000001fe6 <+8166>:	jb     0x2523 <stream_kernel+9507>

711	        norm[5]*norm[22] >= tolerance &&
   0x0000000000001fec <+8172>:	movss  0x2c(%rdx,%rsi,1),%xmm9
   0x0000000000001ff3 <+8179>:	mulss  0x70(%rdx,%rsi,1),%xmm9
   0x0000000000001ffa <+8186>:	comiss %xmm1,%xmm9
   0x0000000000001ffe <+8190>:	jb     0x2523 <stream_kernel+9507>

712	        norm[6]*norm[26] >= tolerance &&
   0x0000000000002004 <+8196>:	movss  0x30(%rdx,%rsi,1),%xmm9
   0x000000000000200b <+8203>:	mulss  0x80(%rdx,%rsi,1),%xmm9
   0x0000000000002015 <+8213>:	comiss %xmm1,%xmm9
   0x0000000000002019 <+8217>:	jb     0x2523 <stream_kernel+9507>

713	        norm[7]*norm[30] >= tolerance)
   0x000000000000201f <+8223>:	movss  0x34(%rdx,%rsi,1),%xmm9
   0x0000000000002026 <+8230>:	mulss  0x90(%rdx,%rsi,1),%xmm9
   0x0000000000002030 <+8240>:	comiss %xmm1,%xmm9
   0x0000000000002034 <+8244>:	jb     0x2523 <stream_kernel+9507>

714	    {
715	      /* A(2,1)*B(1,3) = C(2,3). */
716	      for (i = 0; i < 4; i++)
717	      {
718	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
719	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
720	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000203a <+8250>:	movaps 0x400(%r8),%xmm10
   0x000000000000205a <+8282>:	movaps 0x440(%r8),%xmm13
   0x0000000000002072 <+8306>:	movaps 0x480(%r8),%xmm8
   0x000000000000207a <+8314>:	mulps  0x80(%rcx),%xmm10
   0x000000000000209a <+8346>:	mulps  0x80(%rcx),%xmm13
   0x00000000000020b1 <+8369>:	mulps  0x80(%rcx),%xmm8
   0x0000000000002109 <+8457>:	movaps 0x4c0(%r8),%xmm13
   0x0000000000002111 <+8465>:	mulps  0x80(%rcx),%xmm13

721
722	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
723	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
724	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002042 <+8258>:	movaps 0x410(%r8),%xmm11
   0x0000000000002062 <+8290>:	movaps 0x450(%r8),%xmm14
   0x0000000000002082 <+8322>:	mulps  0x90(%rcx),%xmm11
   0x00000000000020a2 <+8354>:	mulps  0x90(%rcx),%xmm14
   0x00000000000020b9 <+8377>:	movaps 0x4d0(%r8),%xmm15
   0x00000000000020c1 <+8385>:	addps  %xmm10,%xmm11
   0x00000000000020c5 <+8389>:	mulps  0x90(%rcx),%xmm15
   0x00000000000020d9 <+8409>:	movaps 0x490(%r8),%xmm11
   0x00000000000020e9 <+8425>:	mulps  0x90(%rcx),%xmm11
   0x0000000000002105 <+8453>:	addps  %xmm13,%xmm14
   0x0000000000002141 <+8513>:	addps  %xmm8,%xmm11
   0x000000000000217d <+8573>:	addps  %xmm13,%xmm15

725
726	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
727	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
728	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000204a <+8266>:	movaps 0x420(%r8),%xmm12
   0x000000000000208a <+8330>:	mulps  0xa0(%rcx),%xmm12
   0x00000000000020cd <+8397>:	movaps 0x460(%r8),%xmm10
   0x00000000000020d5 <+8405>:	addps  %xmm11,%xmm12
   0x00000000000020e1 <+8417>:	mulps  0xa0(%rcx),%xmm10
   0x00000000000020f5 <+8437>:	movaps 0x4a0(%r8),%xmm12
   0x00000000000020fd <+8445>:	mulps  0xa0(%rcx),%xmm12
   0x0000000000002119 <+8473>:	addps  %xmm14,%xmm10
   0x0000000000002131 <+8497>:	movaps 0x4e0(%r8),%xmm10
   0x0000000000002139 <+8505>:	mulps  0xa0(%rcx),%xmm10
   0x0000000000002155 <+8533>:	addps  %xmm11,%xmm12
   0x0000000000002191 <+8593>:	addps  %xmm15,%xmm10

729
730	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
731	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
732	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002052 <+8274>:	movaps 0x430(%r8),%xmm9
   0x000000000000206a <+8298>:	movaps 0x470(%r8),%xmm7
   0x0000000000002092 <+8338>:	mulps  0xb0(%rcx),%xmm9
   0x00000000000020aa <+8362>:	mulps  0xb0(%rcx),%xmm7
   0x00000000000020f1 <+8433>:	addps  %xmm12,%xmm9
   0x000000000000212d <+8493>:	addps  %xmm10,%xmm7
   0x0000000000002145 <+8517>:	movaps 0x4b0(%r8),%xmm8
   0x000000000000214d <+8525>:	mulps  0xb0(%rcx),%xmm8
   0x0000000000002169 <+8553>:	addps  %xmm12,%xmm8
   0x0000000000002195 <+8597>:	movaps 0x4f0(%r8),%xmm15
   0x000000000000219d <+8605>:	mulps  0xb0(%rcx),%xmm15
   0x00000000000021b9 <+8633>:	addps  %xmm10,%xmm15

733	      }
734
735	      /* A(2,2)*B(2,3) = C(2,3). */
736	      for (i = 0; i < 4; i++)
737	      {
738	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
739	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
740	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002159 <+8537>:	movaps 0x500(%r8),%xmm11
   0x0000000000002161 <+8545>:	mulps  0x180(%rcx),%xmm11
   0x00000000000021a5 <+8613>:	addps  %xmm9,%xmm11
   0x00000000000021d1 <+8657>:	movaps 0x540(%r8),%xmm11
   0x00000000000021d9 <+8665>:	mulps  0x180(%rcx),%xmm11
   0x00000000000021e5 <+8677>:	movaps 0x580(%r8),%xmm9
   0x00000000000021ed <+8685>:	mulps  0x180(%rcx),%xmm9
   0x0000000000002209 <+8713>:	addps  %xmm7,%xmm11
   0x000000000000221c <+8732>:	addps  %xmm8,%xmm9
   0x000000000000225b <+8795>:	movaps 0x5c0(%r8),%xmm10
   0x0000000000002263 <+8803>:	mulps  0x180(%rcx),%xmm10
   0x00000000000022a7 <+8871>:	addps  %xmm15,%xmm10

741
742	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
743	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
744	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000021a9 <+8617>:	movaps 0x510(%r8),%xmm9
   0x00000000000021b1 <+8625>:	mulps  0x190(%rcx),%xmm9
   0x00000000000021cd <+8653>:	addps  %xmm11,%xmm9
   0x000000000000220d <+8717>:	movaps 0x550(%r8),%xmm7
   0x0000000000002215 <+8725>:	mulps  0x190(%rcx),%xmm7
   0x0000000000002220 <+8736>:	movaps 0x590(%r8),%xmm8
   0x0000000000002228 <+8744>:	mulps  0x190(%rcx),%xmm8
   0x0000000000002230 <+8752>:	addps  %xmm11,%xmm7
   0x000000000000226b <+8811>:	addps  %xmm9,%xmm8
   0x00000000000022ab <+8875>:	movaps 0x5d0(%r8),%xmm15
   0x00000000000022b3 <+8883>:	mulps  0x190(%rcx),%xmm15
   0x00000000000022cf <+8911>:	addps  %xmm10,%xmm15

745
746	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
747	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
748	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000021bd <+8637>:	movaps 0x520(%r8),%xmm10
   0x00000000000021c5 <+8645>:	mulps  0x1a0(%rcx),%xmm10
   0x00000000000021e1 <+8673>:	addps  %xmm9,%xmm10
   0x00000000000021f9 <+8697>:	movaps 0x560(%r8),%xmm10
   0x0000000000002201 <+8705>:	mulps  0x1a0(%rcx),%xmm10
   0x0000000000002234 <+8756>:	movaps 0x5a0(%r8),%xmm11
   0x000000000000223c <+8764>:	mulps  0x1a0(%rcx),%xmm11
   0x0000000000002244 <+8772>:	addps  %xmm7,%xmm10
   0x000000000000227f <+8831>:	addps  %xmm8,%xmm11
   0x00000000000022d3 <+8915>:	movaps 0x5e0(%r8),%xmm10
   0x00000000000022db <+8923>:	mulps  0x1a0(%rcx),%xmm10
   0x00000000000022f7 <+8951>:	addps  %xmm15,%xmm10

749
750	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
751	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
752	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000211d <+8477>:	movaps 0x530(%r8),%xmm14
   0x0000000000002125 <+8485>:	mulps  0x1b0(%rcx),%xmm14
   0x000000000000216d <+8557>:	movaps 0x5b0(%r8),%xmm12
   0x0000000000002175 <+8565>:	mulps  0x1b0(%rcx),%xmm12
   0x0000000000002181 <+8577>:	movaps 0x570(%r8),%xmm13
   0x0000000000002189 <+8585>:	mulps  0x1b0(%rcx),%xmm13
   0x00000000000021f5 <+8693>:	addps  %xmm10,%xmm14
   0x0000000000002257 <+8791>:	addps  %xmm10,%xmm13
   0x0000000000002293 <+8851>:	addps  %xmm11,%xmm12
   0x0000000000002297 <+8855>:	movaps 0x5f0(%r8),%xmm11
   0x000000000000229f <+8863>:	mulps  0x1b0(%rcx),%xmm11
   0x000000000000230b <+8971>:	addps  %xmm10,%xmm11

753	      }
754
755	      /* A(2,3)*B(3,3) = C(2,3). */
756	      for (i = 0; i < 4; i++)
757	      {
758	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
759	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
760	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002248 <+8776>:	movaps 0x640(%r8),%xmm7
   0x0000000000002250 <+8784>:	mulps  0x280(%rcx),%xmm7
   0x000000000000226f <+8815>:	movaps 0x680(%r8),%xmm9
   0x0000000000002277 <+8823>:	mulps  0x280(%rcx),%xmm9
   0x0000000000002283 <+8835>:	movaps 0x600(%r8),%xmm8
   0x000000000000228b <+8843>:	mulps  0x280(%rcx),%xmm8
   0x00000000000022bb <+8891>:	addps  %xmm14,%xmm8
   0x0000000000002347 <+9031>:	addps  %xmm13,%xmm7
   0x000000000000235b <+9051>:	addps  %xmm12,%xmm9
   0x000000000000239a <+9114>:	movaps 0x6c0(%r8),%xmm13
   0x00000000000023a2 <+9122>:	mulps  0x280(%rcx),%xmm13
   0x00000000000023e5 <+9189>:	addps  %xmm11,%xmm13

761
762	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
763	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
764	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000022bf <+8895>:	movaps 0x610(%r8),%xmm14
   0x00000000000022c7 <+8903>:	mulps  0x290(%rcx),%xmm14
   0x00000000000022e3 <+8931>:	addps  %xmm8,%xmm14
   0x000000000000234b <+9035>:	movaps 0x650(%r8),%xmm13
   0x0000000000002353 <+9043>:	mulps  0x290(%rcx),%xmm13
   0x000000000000235f <+9055>:	movaps 0x690(%r8),%xmm12
   0x0000000000002367 <+9063>:	mulps  0x290(%rcx),%xmm12
   0x000000000000236f <+9071>:	addps  %xmm7,%xmm13
   0x0000000000002382 <+9090>:	addps  %xmm9,%xmm12
   0x00000000000023e9 <+9193>:	movaps 0x6d0(%r8),%xmm11
   0x00000000000023f1 <+9201>:	mulps  0x290(%rcx),%xmm11
   0x000000000000240d <+9229>:	addps  %xmm13,%xmm11

765
766	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
767	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
768	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000022e7 <+8935>:	movaps 0x620(%r8),%xmm8
   0x00000000000022ef <+8943>:	mulps  0x2a0(%rcx),%xmm8
   0x000000000000231f <+8991>:	addps  %xmm14,%xmm8
   0x0000000000002323 <+8995>:	movaps 0x6e0(%r8),%xmm14
   0x000000000000232b <+9003>:	mulps  0x2a0(%rcx),%xmm14
   0x0000000000002373 <+9075>:	movaps 0x660(%r8),%xmm7
   0x000000000000237b <+9083>:	mulps  0x2a0(%rcx),%xmm7
   0x0000000000002386 <+9094>:	movaps 0x6a0(%r8),%xmm9
   0x000000000000238e <+9102>:	mulps  0x2a0(%rcx),%xmm9
   0x0000000000002396 <+9110>:	addps  %xmm13,%xmm7
   0x00000000000023bd <+9149>:	addps  %xmm12,%xmm9
   0x0000000000002421 <+9249>:	addps  %xmm11,%xmm14

769
770	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
771	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
772	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000230f <+8975>:	movaps 0x630(%r8),%xmm10
   0x0000000000002317 <+8983>:	mulps  0x2b0(%rcx),%xmm10
   0x0000000000002333 <+9011>:	addps  %xmm8,%xmm10
   0x0000000000002337 <+9015>:	movaps 0x670(%r8),%xmm8
   0x000000000000233f <+9023>:	mulps  0x2b0(%rcx),%xmm8
   0x00000000000023aa <+9130>:	addps  %xmm7,%xmm8
   0x00000000000023ae <+9134>:	movaps 0x6b0(%r8),%xmm7
   0x00000000000023b6 <+9142>:	mulps  0x2b0(%rcx),%xmm7
   0x00000000000023d1 <+9169>:	addps  %xmm9,%xmm7
   0x00000000000023d5 <+9173>:	movaps 0x6f0(%r8),%xmm9
   0x00000000000023dd <+9181>:	mulps  0x2b0(%rcx),%xmm9
   0x0000000000002435 <+9269>:	addps  %xmm14,%xmm9

773	      }
774
775	      /* A(2,4)*B(4,3) = C(2,3). */
776	      for (i = 0; i < 4; i++)
777	      {
778	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
779	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
780	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000022fb <+8955>:	movaps 0x700(%r8),%xmm15
   0x0000000000002303 <+8963>:	mulps  0x380(%rcx),%xmm15
   0x00000000000023c1 <+9153>:	movaps 0x740(%r8),%xmm12
   0x00000000000023c9 <+9161>:	mulps  0x380(%rcx),%xmm12
   0x00000000000023f9 <+9209>:	addps  %xmm10,%xmm15
   0x0000000000002439 <+9273>:	movaps 0x780(%r8),%xmm14
   0x0000000000002441 <+9281>:	mulps  0x380(%rcx),%xmm14
   0x0000000000002461 <+9313>:	addps  %xmm8,%xmm12
   0x000000000000248d <+9357>:	movaps 0x7c0(%r8),%xmm12
   0x0000000000002495 <+9365>:	mulps  0x380(%rcx),%xmm12
   0x00000000000024b1 <+9393>:	addps  %xmm7,%xmm14
   0x00000000000024c8 <+9416>:	addps  %xmm9,%xmm12

781
782	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
783	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
784	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000023fd <+9213>:	movaps 0x710(%r8),%xmm10
   0x0000000000002405 <+9221>:	mulps  0x390(%rcx),%xmm10
   0x0000000000002449 <+9289>:	addps  %xmm15,%xmm10
   0x0000000000002465 <+9317>:	movaps 0x750(%r8),%xmm8
   0x000000000000246d <+9325>:	mulps  0x390(%rcx),%xmm8
   0x0000000000002489 <+9353>:	addps  %xmm12,%xmm8
   0x00000000000024b5 <+9397>:	movaps 0x790(%r8),%xmm7
   0x00000000000024bd <+9405>:	mulps  0x390(%rcx),%xmm7
   0x00000000000024cc <+9420>:	movaps 0x7d0(%r8),%xmm9
   0x00000000000024d4 <+9428>:	mulps  0x390(%rcx),%xmm9
   0x00000000000024dc <+9436>:	addps  %xmm14,%xmm7
   0x00000000000024f3 <+9459>:	addps  %xmm12,%xmm9

785
786	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
787	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
788	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002411 <+9233>:	movaps 0x760(%r8),%xmm13
   0x0000000000002419 <+9241>:	mulps  0x3a0(%rcx),%xmm13
   0x0000000000002425 <+9253>:	movaps 0x720(%r8),%xmm11
   0x000000000000242d <+9261>:	mulps  0x3a0(%rcx),%xmm11
   0x000000000000244d <+9293>:	addps  %xmm10,%xmm11
   0x0000000000002479 <+9337>:	movaps 0x7a0(%r8),%xmm11
   0x0000000000002481 <+9345>:	mulps  0x3a0(%rcx),%xmm11
   0x000000000000249d <+9373>:	addps  %xmm8,%xmm13
   0x00000000000024e0 <+9440>:	addps  %xmm7,%xmm11
   0x000000000000250b <+9483>:	movaps 0x7e0(%r8),%xmm11
   0x0000000000002513 <+9491>:	mulps  0x3a0(%rcx),%xmm11
   0x000000000000251b <+9499>:	addps  %xmm9,%xmm11

789
790	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
791	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
792	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002451 <+9297>:	movaps 0x730(%r8),%xmm10
   0x0000000000002459 <+9305>:	mulps  0x3b0(%rcx),%xmm10
   0x0000000000002475 <+9333>:	addps  %xmm11,%xmm10
   0x00000000000024a1 <+9377>:	movaps 0x770(%r8),%xmm8
   0x00000000000024a9 <+9385>:	mulps  0x3b0(%rcx),%xmm8
   0x00000000000024c4 <+9412>:	addps  %xmm13,%xmm8
   0x00000000000024e4 <+9444>:	movaps 0x7b0(%r8),%xmm7
   0x00000000000024ec <+9452>:	mulps  0x3b0(%rcx),%xmm7
   0x00000000000024f7 <+9463>:	movaps 0x7f0(%r8),%xmm12
   0x00000000000024ff <+9471>:	mulps  0x3b0(%rcx),%xmm12
   0x0000000000002507 <+9479>:	addps  %xmm11,%xmm7
   0x000000000000251f <+9503>:	addps  %xmm11,%xmm12

793	      }
794	    }
795
796	    /* Store C(2,3) block. */
797	    for (i = 0; i < 4; i++)
798	    {
799	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000002523 <+9507>:	mulps  %xmm6,%xmm10
   0x000000000000253f <+9535>:	mulps  %xmm6,%xmm8
   0x0000000000002557 <+9559>:	mulps  %xmm6,%xmm7
   0x000000000000256d <+9581>:	mulps  %xmm6,%xmm12

800	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_23]), C_row[i]);
   0x0000000000002527 <+9511>:	addps  0x180(%r9),%xmm10
   0x0000000000002543 <+9539>:	addps  0x190(%r9),%xmm8
   0x000000000000255a <+9562>:	addps  0x1a0(%r9),%xmm7
   0x0000000000002571 <+9585>:	addps  0x1b0(%r9),%xmm12

801	      _mm_store_ps(&C[i*4+C_OFFSET_23], C_row[i]);
   0x0000000000002537 <+9527>:	movaps %xmm10,0x180(%r9)
   0x000000000000254b <+9547>:	movaps %xmm8,0x190(%r9)
   0x0000000000002562 <+9570>:	movaps %xmm7,0x1a0(%r9)
   0x0000000000002579 <+9593>:	movaps %xmm12,0x1b0(%r9)

802	    }
803
804	    /* Reset C(2,4) matrix accumulators */
805	    C_row[0] = _mm_setzero_ps();
   0x000000000000252f <+9519>:	xorps  %xmm9,%xmm9

806	    C_row[1] = _mm_setzero_ps();
   0x0000000000002553 <+9555>:	xorps  %xmm8,%xmm8

807	    C_row[2] = _mm_setzero_ps();
   0x000000000000256a <+9578>:	xorps  %xmm7,%xmm7

808	    C_row[3] = _mm_setzero_ps();
   0x0000000000002533 <+9523>:	xorps  %xmm11,%xmm11

809
810	    if (norm[4]*norm[19] >= tolerance &&
   0x0000000000002581 <+9601>:	mulss  %xmm2,%xmm0
   0x0000000000002585 <+9605>:	comiss %xmm1,%xmm0
   0x0000000000002588 <+9608>:	jb     0x2aad <stream_kernel+10925>

811	        norm[5]*norm[23] >= tolerance &&
   0x000000000000258e <+9614>:	movss  0x2c(%rdx,%rsi,1),%xmm0
   0x0000000000002594 <+9620>:	mulss  0x74(%rdx,%rsi,1),%xmm0
   0x000000000000259a <+9626>:	comiss %xmm1,%xmm0
   0x000000000000259d <+9629>:	jb     0x2aad <stream_kernel+10925>

812	        norm[6]*norm[27] >= tolerance &&
   0x00000000000025a3 <+9635>:	movss  0x30(%rdx,%rsi,1),%xmm0
   0x00000000000025a9 <+9641>:	mulss  0x84(%rdx,%rsi,1),%xmm0
   0x00000000000025b2 <+9650>:	comiss %xmm1,%xmm0
   0x00000000000025b5 <+9653>:	jb     0x2aad <stream_kernel+10925>

813	        norm[7]*norm[31] >= tolerance)
   0x00000000000025bb <+9659>:	movss  0x34(%rdx,%rsi,1),%xmm0
   0x00000000000025c1 <+9665>:	mulss  0x94(%rdx,%rsi,1),%xmm0
   0x00000000000025ca <+9674>:	comiss %xmm1,%xmm0
   0x00000000000025cd <+9677>:	jb     0x2aad <stream_kernel+10925>

814	    {
815	      /* A(2,1)*B(1,4) = C(2,4). */
816	      for (i = 0; i < 4; i++)
817	      {
818	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
819	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
820	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000025d3 <+9683>:	movaps 0x400(%r8),%xmm7
   0x00000000000025f3 <+9715>:	movaps 0x440(%r8),%xmm11
   0x0000000000002613 <+9747>:	movaps 0x480(%r8),%xmm14
   0x0000000000002623 <+9763>:	mulps  0xc0(%rcx),%xmm7
   0x0000000000002641 <+9793>:	mulps  0xc0(%rcx),%xmm11
   0x0000000000002661 <+9825>:	mulps  0xc0(%rcx),%xmm14
   0x000000000000269a <+9882>:	movaps 0x4c0(%r8),%xmm10
   0x00000000000026a2 <+9890>:	mulps  0xc0(%rcx),%xmm10

821
822	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
823	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
824	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000025db <+9691>:	movaps 0x410(%r8),%xmm0
   0x00000000000025fb <+9723>:	movaps 0x450(%r8),%xmm12
   0x000000000000262a <+9770>:	mulps  0xd0(%rcx),%xmm0
   0x0000000000002649 <+9801>:	mulps  0xd0(%rcx),%xmm12
   0x0000000000002671 <+9841>:	addps  %xmm7,%xmm0
   0x0000000000002674 <+9844>:	movaps 0x490(%r8),%xmm7
   0x000000000000267c <+9852>:	mulps  0xd0(%rcx),%xmm7
   0x00000000000026aa <+9898>:	addps  %xmm11,%xmm12
   0x00000000000026ae <+9902>:	movaps 0x4d0(%r8),%xmm11
   0x00000000000026b6 <+9910>:	mulps  0xd0(%rcx),%xmm11
   0x00000000000026e6 <+9958>:	addps  %xmm14,%xmm7
   0x000000000000271f <+10015>:	addps  %xmm10,%xmm11

825
826	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
827	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
828	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000025e3 <+9699>:	movaps 0x420(%r8),%xmm10
   0x0000000000002603 <+9731>:	movaps 0x460(%r8),%xmm13
   0x0000000000002631 <+9777>:	mulps  0xe0(%rcx),%xmm10
   0x0000000000002651 <+9809>:	mulps  0xe0(%rcx),%xmm13
   0x0000000000002683 <+9859>:	addps  %xmm0,%xmm10
   0x0000000000002687 <+9863>:	movaps 0x4a0(%r8),%xmm0
   0x000000000000268f <+9871>:	mulps  0xe0(%rcx),%xmm0
   0x00000000000026be <+9918>:	addps  %xmm12,%xmm13
   0x00000000000026c2 <+9922>:	movaps 0x4e0(%r8),%xmm12
   0x00000000000026ca <+9930>:	mulps  0xe0(%rcx),%xmm12
   0x00000000000026fa <+9978>:	addps  %xmm7,%xmm0
   0x0000000000002733 <+10035>:	addps  %xmm11,%xmm12

829
830	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
831	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
832	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000025eb <+9707>:	movaps 0x430(%r8),%xmm8
   0x000000000000260b <+9739>:	movaps 0x470(%r8),%xmm9
   0x000000000000261b <+9755>:	movaps 0x4b0(%r8),%xmm15
   0x0000000000002639 <+9785>:	mulps  0xf0(%rcx),%xmm8
   0x0000000000002659 <+9817>:	mulps  0xf0(%rcx),%xmm9
   0x0000000000002669 <+9833>:	mulps  0xf0(%rcx),%xmm15
   0x0000000000002696 <+9878>:	addps  %xmm10,%xmm8
   0x00000000000026d2 <+9938>:	addps  %xmm13,%xmm9
   0x00000000000026ea <+9962>:	movaps 0x4f0(%r8),%xmm14
   0x00000000000026f2 <+9970>:	mulps  0xf0(%rcx),%xmm14
   0x000000000000270c <+9996>:	addps  %xmm0,%xmm15
   0x0000000000002747 <+10055>:	addps  %xmm12,%xmm14

833	      }
834
835	      /* A(2,2)*B(2,4) = C(2,4). */
836	      for (i = 0; i < 4; i++)
837	      {
838	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
839	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
840	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000026fd <+9981>:	movaps 0x500(%r8),%xmm7
   0x0000000000002705 <+9989>:	mulps  0x1c0(%rcx),%xmm7
   0x0000000000002710 <+10000>:	movaps 0x540(%r8),%xmm0
   0x0000000000002718 <+10008>:	mulps  0x1c0(%rcx),%xmm0
   0x000000000000275b <+10075>:	addps  %xmm8,%xmm7
   0x000000000000276f <+10095>:	addps  %xmm9,%xmm0
   0x00000000000027c1 <+10177>:	movaps 0x580(%r8),%xmm7
   0x00000000000027c9 <+10185>:	mulps  0x1c0(%rcx),%xmm7
   0x00000000000027e8 <+10216>:	movaps 0x5c0(%r8),%xmm8
   0x00000000000027f0 <+10224>:	mulps  0x1c0(%rcx),%xmm8
   0x00000000000027f8 <+10232>:	addps  %xmm15,%xmm7
   0x000000000000280c <+10252>:	addps  %xmm14,%xmm8

841
842	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
843	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
844	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000275f <+10079>:	movaps 0x510(%r8),%xmm8
   0x0000000000002767 <+10087>:	mulps  0x1d0(%rcx),%xmm8
   0x0000000000002773 <+10099>:	movaps 0x550(%r8),%xmm9
   0x000000000000277b <+10107>:	mulps  0x1d0(%rcx),%xmm9
   0x0000000000002783 <+10115>:	addps  %xmm7,%xmm8
   0x0000000000002796 <+10134>:	addps  %xmm0,%xmm9
   0x00000000000027fc <+10236>:	movaps 0x590(%r8),%xmm15
   0x0000000000002804 <+10244>:	mulps  0x1d0(%rcx),%xmm15
   0x0000000000002820 <+10272>:	addps  %xmm7,%xmm15
   0x0000000000002824 <+10276>:	movaps 0x5d0(%r8),%xmm7
   0x000000000000282c <+10284>:	mulps  0x1d0(%rcx),%xmm7
   0x000000000000285b <+10331>:	addps  %xmm8,%xmm7

845
846	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
847	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
848	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002787 <+10119>:	movaps 0x520(%r8),%xmm7
   0x000000000000278f <+10127>:	mulps  0x1e0(%rcx),%xmm7
   0x000000000000279a <+10138>:	movaps 0x5e0(%r8),%xmm0
   0x00000000000027a2 <+10146>:	mulps  0x1e0(%rcx),%xmm0
   0x00000000000027a9 <+10153>:	addps  %xmm8,%xmm7
   0x00000000000027ad <+10157>:	movaps 0x560(%r8),%xmm8
   0x00000000000027b5 <+10165>:	mulps  0x1e0(%rcx),%xmm8
   0x00000000000027d0 <+10192>:	addps  %xmm9,%xmm8
   0x00000000000027d4 <+10196>:	movaps 0x5a0(%r8),%xmm9
   0x00000000000027dc <+10204>:	mulps  0x1e0(%rcx),%xmm9
   0x0000000000002833 <+10291>:	addps  %xmm15,%xmm9
   0x000000000000286f <+10351>:	addps  %xmm7,%xmm0

849
850	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
851	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
852	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000026d6 <+9942>:	movaps 0x530(%r8),%xmm13
   0x00000000000026de <+9950>:	mulps  0x1f0(%rcx),%xmm13
   0x0000000000002723 <+10019>:	movaps 0x5f0(%r8),%xmm10
   0x000000000000272b <+10027>:	mulps  0x1f0(%rcx),%xmm10
   0x0000000000002737 <+10039>:	movaps 0x5b0(%r8),%xmm11
   0x000000000000273f <+10047>:	mulps  0x1f0(%rcx),%xmm11
   0x000000000000274b <+10059>:	movaps 0x570(%r8),%xmm12
   0x0000000000002753 <+10067>:	mulps  0x1f0(%rcx),%xmm12
   0x00000000000027bd <+10173>:	addps  %xmm7,%xmm13
   0x00000000000027e4 <+10212>:	addps  %xmm8,%xmm12
   0x0000000000002847 <+10311>:	addps  %xmm9,%xmm11
   0x0000000000002881 <+10369>:	addps  %xmm0,%xmm10

853	      }
854
855	      /* A(2,3)*B(3,4) = C(2,4). */
856	      for (i = 0; i < 4; i++)
857	      {
858	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_23]);
859	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
860	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002810 <+10256>:	movaps 0x600(%r8),%xmm14
   0x0000000000002818 <+10264>:	mulps  0x2c0(%rcx),%xmm14
   0x0000000000002872 <+10354>:	movaps 0x640(%r8),%xmm7
   0x000000000000287a <+10362>:	mulps  0x2c0(%rcx),%xmm7
   0x0000000000002894 <+10388>:	addps  %xmm13,%xmm14
   0x00000000000028a8 <+10408>:	addps  %xmm12,%xmm7
   0x00000000000028c0 <+10432>:	movaps 0x680(%r8),%xmm14
   0x00000000000028c8 <+10440>:	mulps  0x2c0(%rcx),%xmm14
   0x000000000000290e <+10510>:	movaps 0x6c0(%r8),%xmm0
   0x0000000000002916 <+10518>:	mulps  0x2c0(%rcx),%xmm0
   0x0000000000002931 <+10545>:	addps  %xmm11,%xmm14
   0x0000000000002945 <+10565>:	addps  %xmm10,%xmm0

861
862	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_23]);
863	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
864	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002885 <+10373>:	movaps 0x650(%r8),%xmm0
   0x000000000000288d <+10381>:	mulps  0x2d0(%rcx),%xmm0
   0x0000000000002898 <+10392>:	movaps 0x610(%r8),%xmm13
   0x00000000000028a0 <+10400>:	mulps  0x2d0(%rcx),%xmm13
   0x00000000000028bc <+10428>:	addps  %xmm14,%xmm13
   0x00000000000028f8 <+10488>:	addps  %xmm7,%xmm0
   0x0000000000002935 <+10549>:	movaps 0x690(%r8),%xmm11
   0x000000000000293d <+10557>:	mulps  0x2d0(%rcx),%xmm11
   0x0000000000002949 <+10569>:	movaps 0x6d0(%r8),%xmm10
   0x0000000000002951 <+10577>:	mulps  0x2d0(%rcx),%xmm10
   0x0000000000002959 <+10585>:	addps  %xmm14,%xmm11
   0x0000000000002995 <+10645>:	addps  %xmm0,%xmm10

865
866	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_23]);
867	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
868	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002837 <+10295>:	movaps 0x6a0(%r8),%xmm15
   0x000000000000283f <+10303>:	mulps  0x2e0(%rcx),%xmm15
   0x000000000000285f <+10335>:	movaps 0x620(%r8),%xmm8
   0x0000000000002867 <+10343>:	mulps  0x2e0(%rcx),%xmm8
   0x00000000000028ac <+10412>:	movaps 0x660(%r8),%xmm12
   0x00000000000028b4 <+10420>:	mulps  0x2e0(%rcx),%xmm12
   0x00000000000028d0 <+10448>:	addps  %xmm13,%xmm8
   0x000000000000290a <+10506>:	addps  %xmm0,%xmm12
   0x000000000000296d <+10605>:	addps  %xmm11,%xmm15
   0x0000000000002971 <+10609>:	movaps 0x6e0(%r8),%xmm11
   0x0000000000002979 <+10617>:	mulps  0x2e0(%rcx),%xmm11
   0x00000000000029a8 <+10664>:	addps  %xmm10,%xmm11

869
870	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_23]);
871	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
872	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000284b <+10315>:	movaps 0x630(%r8),%xmm9
   0x0000000000002853 <+10323>:	mulps  0x2f0(%rcx),%xmm9
   0x00000000000028e4 <+10468>:	addps  %xmm8,%xmm9
   0x00000000000028e8 <+10472>:	movaps 0x670(%r8),%xmm8
   0x00000000000028f0 <+10480>:	mulps  0x2f0(%rcx),%xmm8
   0x00000000000028fb <+10491>:	movaps 0x6b0(%r8),%xmm7
   0x0000000000002903 <+10499>:	mulps  0x2f0(%rcx),%xmm7
   0x000000000000291d <+10525>:	addps  %xmm12,%xmm8
   0x0000000000002981 <+10625>:	addps  %xmm15,%xmm7
   0x0000000000002999 <+10649>:	movaps 0x6f0(%r8),%xmm0
   0x00000000000029a1 <+10657>:	mulps  0x2f0(%rcx),%xmm0
   0x00000000000029bc <+10684>:	addps  %xmm11,%xmm0

873	      }
874
875	      /* A(2,4)*B(4,4) = C(2,4). */
876	      for (i = 0; i < 4; i++)
877	      {
878	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_24]);
879	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
880	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002921 <+10529>:	movaps 0x700(%r8),%xmm12
   0x0000000000002929 <+10537>:	mulps  0x3c0(%rcx),%xmm12
   0x000000000000295d <+10589>:	movaps 0x740(%r8),%xmm14
   0x0000000000002965 <+10597>:	mulps  0x3c0(%rcx),%xmm14
   0x00000000000029ac <+10668>:	movaps 0x780(%r8),%xmm10
   0x00000000000029b4 <+10676>:	mulps  0x3c0(%rcx),%xmm10
   0x00000000000029d0 <+10704>:	addps  %xmm9,%xmm12
   0x00000000000029e4 <+10724>:	addps  %xmm8,%xmm14
   0x00000000000029fc <+10748>:	movaps 0x7c0(%r8),%xmm12
   0x0000000000002a04 <+10756>:	mulps  0x3c0(%rcx),%xmm12
   0x0000000000002a3c <+10812>:	addps  %xmm7,%xmm10
   0x0000000000002a53 <+10835>:	addps  %xmm0,%xmm12

881
882	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_24]);
883	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
884	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000029d4 <+10708>:	movaps 0x710(%r8),%xmm9
   0x00000000000029dc <+10716>:	mulps  0x3d0(%rcx),%xmm9
   0x00000000000029e8 <+10728>:	movaps 0x750(%r8),%xmm8
   0x00000000000029f0 <+10736>:	mulps  0x3d0(%rcx),%xmm8
   0x00000000000029f8 <+10744>:	addps  %xmm12,%xmm9
   0x0000000000002a20 <+10784>:	addps  %xmm14,%xmm8
   0x0000000000002a40 <+10816>:	movaps 0x790(%r8),%xmm7
   0x0000000000002a48 <+10824>:	mulps  0x3d0(%rcx),%xmm7
   0x0000000000002a57 <+10839>:	movaps 0x7d0(%r8),%xmm0
   0x0000000000002a5f <+10847>:	mulps  0x3d0(%rcx),%xmm0
   0x0000000000002a66 <+10854>:	addps  %xmm10,%xmm7
   0x0000000000002a8d <+10893>:	addps  %xmm12,%xmm0

885
886	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_24]);
887	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
888	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000028d4 <+10452>:	movaps 0x720(%r8),%xmm13
   0x00000000000028dc <+10460>:	mulps  0x3e0(%rcx),%xmm13
   0x0000000000002985 <+10629>:	movaps 0x760(%r8),%xmm15
   0x000000000000298d <+10637>:	mulps  0x3e0(%rcx),%xmm15
   0x00000000000029c0 <+10688>:	movaps 0x7a0(%r8),%xmm11
   0x00000000000029c8 <+10696>:	mulps  0x3e0(%rcx),%xmm11
   0x0000000000002a0c <+10764>:	addps  %xmm9,%xmm13
   0x0000000000002a24 <+10788>:	addps  %xmm8,%xmm15
   0x0000000000002a6a <+10858>:	movaps 0x7e0(%r8),%xmm10
   0x0000000000002a72 <+10866>:	mulps  0x3e0(%rcx),%xmm10
   0x0000000000002a7a <+10874>:	addps  %xmm7,%xmm11
   0x0000000000002a91 <+10897>:	addps  %xmm0,%xmm10

889
890	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_24]);
891	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
892	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002a10 <+10768>:	movaps 0x730(%r8),%xmm9
   0x0000000000002a18 <+10776>:	mulps  0x3f0(%rcx),%xmm9
   0x0000000000002a28 <+10792>:	movaps 0x770(%r8),%xmm8
   0x0000000000002a30 <+10800>:	mulps  0x3f0(%rcx),%xmm8
   0x0000000000002a38 <+10808>:	addps  %xmm13,%xmm9
   0x0000000000002a4f <+10831>:	addps  %xmm15,%xmm8
   0x0000000000002a7e <+10878>:	movaps 0x7b0(%r8),%xmm7
   0x0000000000002a86 <+10886>:	mulps  0x3f0(%rcx),%xmm7
   0x0000000000002a95 <+10901>:	addps  %xmm11,%xmm7
   0x0000000000002a99 <+10905>:	movaps 0x7f0(%r8),%xmm11
   0x0000000000002aa1 <+10913>:	mulps  0x3f0(%rcx),%xmm11
   0x0000000000002aa9 <+10921>:	addps  %xmm10,%xmm11

893	      }
894	    }
895
896	    /* Store C(2,4) block. */
897	    for (i = 0; i < 4; i++)
898	    {
899	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000002aad <+10925>:	mulps  %xmm6,%xmm9
   0x0000000000002ad3 <+10963>:	mulps  %xmm6,%xmm8
   0x0000000000002aeb <+10987>:	mulps  %xmm6,%xmm7
   0x0000000000002b01 <+11009>:	mulps  %xmm6,%xmm11

900	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_24]), C_row[i]);
   0x0000000000002ab1 <+10929>:	addps  0x1c0(%r9),%xmm9
   0x0000000000002ad7 <+10967>:	addps  0x1d0(%r9),%xmm8
   0x0000000000002aee <+10990>:	addps  0x1e0(%r9),%xmm7
   0x0000000000002b05 <+11013>:	addps  0x1f0(%r9),%xmm11

901	      _mm_store_ps(&C[i*4+C_OFFSET_24], C_row[i]);
   0x0000000000002ac3 <+10947>:	movaps %xmm9,0x1c0(%r9)
   0x0000000000002adf <+10975>:	movaps %xmm8,0x1d0(%r9)
   0x0000000000002af6 <+10998>:	movaps %xmm7,0x1e0(%r9)
   0x0000000000002b0d <+11021>:	movaps %xmm11,0x1f0(%r9)

902	    }
903
904	    /* Reset C(3,1) matrix accumulators */
905	    C_row[0] = _mm_setzero_ps();
   0x0000000000002abf <+10943>:	xorps  %xmm10,%xmm10

906	    C_row[1] = _mm_setzero_ps();
   0x0000000000002ae7 <+10983>:	xorps  %xmm8,%xmm8

907	    C_row[2] = _mm_setzero_ps();
   0x0000000000002afe <+11006>:	xorps  %xmm7,%xmm7

908	    C_row[3] = _mm_setzero_ps();
   0x0000000000002acb <+10955>:	xorps  %xmm12,%xmm12

909
910	    if (norm[8]*norm[16] >= tolerance &&
   0x0000000000002ab9 <+10937>:	movss  0x38(%rdx,%rsi,1),%xmm0
   0x0000000000002acf <+10959>:	movaps %xmm5,%xmm9
   0x0000000000002b15 <+11029>:	mulss  %xmm0,%xmm9
   0x0000000000002b1a <+11034>:	comiss %xmm1,%xmm9
   0x0000000000002b1e <+11038>:	jb     0x3024 <stream_kernel+12324>

911	        norm[9]*norm[20] >= tolerance &&
   0x0000000000002b24 <+11044>:	movss  0x3c(%rdx,%rsi,1),%xmm9
   0x0000000000002b2b <+11051>:	mulss  0x68(%rdx,%rsi,1),%xmm9
   0x0000000000002b32 <+11058>:	comiss %xmm1,%xmm9
   0x0000000000002b36 <+11062>:	jb     0x3024 <stream_kernel+12324>

912	        norm[10]*norm[24] >= tolerance &&
   0x0000000000002b3c <+11068>:	movss  0x40(%rdx,%rsi,1),%xmm9
   0x0000000000002b43 <+11075>:	mulss  0x78(%rdx,%rsi,1),%xmm9
   0x0000000000002b4a <+11082>:	comiss %xmm1,%xmm9
   0x0000000000002b4e <+11086>:	jb     0x3024 <stream_kernel+12324>

913	        norm[11]*norm[28] >= tolerance)
   0x0000000000002b54 <+11092>:	movss  0x44(%rdx,%rsi,1),%xmm9
   0x0000000000002b5b <+11099>:	mulss  0x88(%rdx,%rsi,1),%xmm9
   0x0000000000002b65 <+11109>:	comiss %xmm1,%xmm9
   0x0000000000002b69 <+11113>:	jb     0x3024 <stream_kernel+12324>

914	    {
915	      /* A(3,1)*B(1,1) = C(3,1). */
916	      for (i = 0; i < 4; i++)
917	      {
918	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
919	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
920	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b6f <+11119>:	movaps 0x800(%r8),%xmm10
   0x0000000000002b8f <+11151>:	movaps 0x840(%r8),%xmm13
   0x0000000000002ba7 <+11175>:	movaps 0x880(%r8),%xmm8
   0x0000000000002baf <+11183>:	mulps  (%rcx),%xmm10
   0x0000000000002bc2 <+11202>:	mulps  (%rcx),%xmm13
   0x0000000000002bcf <+11215>:	mulps  (%rcx),%xmm8
   0x0000000000002c17 <+11287>:	movaps 0x8c0(%r8),%xmm13
   0x0000000000002c1f <+11295>:	mulps  (%rcx),%xmm13

921
922	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
923	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
924	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b77 <+11127>:	movaps 0x810(%r8),%xmm11
   0x0000000000002b97 <+11159>:	movaps 0x850(%r8),%xmm14
   0x0000000000002bb3 <+11187>:	mulps  0x10(%rcx),%xmm11
   0x0000000000002bc6 <+11206>:	mulps  0x10(%rcx),%xmm14
   0x0000000000002bd3 <+11219>:	movaps 0x8d0(%r8),%xmm15
   0x0000000000002bdb <+11227>:	addps  %xmm10,%xmm11
   0x0000000000002bdf <+11231>:	mulps  0x10(%rcx),%xmm15
   0x0000000000002bf0 <+11248>:	movaps 0x890(%r8),%xmm11
   0x0000000000002bfd <+11261>:	mulps  0x10(%rcx),%xmm11
   0x0000000000002c13 <+11283>:	addps  %xmm13,%xmm14
   0x0000000000002c48 <+11336>:	addps  %xmm8,%xmm11
   0x0000000000002c81 <+11393>:	addps  %xmm13,%xmm15

925
926	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
927	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
928	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b7f <+11135>:	movaps 0x820(%r8),%xmm12
   0x0000000000002bb8 <+11192>:	mulps  0x20(%rcx),%xmm12
   0x0000000000002be4 <+11236>:	movaps 0x860(%r8),%xmm10
   0x0000000000002bec <+11244>:	addps  %xmm11,%xmm12
   0x0000000000002bf8 <+11256>:	mulps  0x20(%rcx),%xmm10
   0x0000000000002c06 <+11270>:	movaps 0x8a0(%r8),%xmm12
   0x0000000000002c0e <+11278>:	mulps  0x20(%rcx),%xmm12
   0x0000000000002c23 <+11299>:	addps  %xmm14,%xmm10
   0x0000000000002c3b <+11323>:	movaps 0x8e0(%r8),%xmm10
   0x0000000000002c43 <+11331>:	mulps  0x20(%rcx),%xmm10
   0x0000000000002c59 <+11353>:	addps  %xmm11,%xmm12
   0x0000000000002c95 <+11413>:	addps  %xmm15,%xmm10

929
930	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
931	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
932	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002b87 <+11143>:	movaps 0x830(%r8),%xmm9
   0x0000000000002b9f <+11167>:	movaps 0x870(%r8),%xmm7
   0x0000000000002bbd <+11197>:	mulps  0x30(%rcx),%xmm9
   0x0000000000002bcb <+11211>:	mulps  0x30(%rcx),%xmm7
   0x0000000000002c02 <+11266>:	addps  %xmm12,%xmm9
   0x0000000000002c37 <+11319>:	addps  %xmm10,%xmm7
   0x0000000000002c4c <+11340>:	movaps 0x8b0(%r8),%xmm8
   0x0000000000002c54 <+11348>:	mulps  0x30(%rcx),%xmm8
   0x0000000000002c6d <+11373>:	addps  %xmm12,%xmm8
   0x0000000000002c99 <+11417>:	movaps 0x8f0(%r8),%xmm15
   0x0000000000002ca1 <+11425>:	mulps  0x30(%rcx),%xmm15
   0x0000000000002cba <+11450>:	addps  %xmm10,%xmm15

933	      }
934
935	      /* A(3,2)*B(2,1) = C(3,1). */
936	      for (i = 0; i < 4; i++)
937	      {
938	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
939	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
940	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002c5d <+11357>:	movaps 0x900(%r8),%xmm11
   0x0000000000002c65 <+11365>:	mulps  0x100(%rcx),%xmm11
   0x0000000000002ca6 <+11430>:	addps  %xmm9,%xmm11
   0x0000000000002cd2 <+11474>:	movaps 0x940(%r8),%xmm11
   0x0000000000002cda <+11482>:	mulps  0x100(%rcx),%xmm11
   0x0000000000002ce6 <+11494>:	movaps 0x980(%r8),%xmm9
   0x0000000000002cee <+11502>:	mulps  0x100(%rcx),%xmm9
   0x0000000000002d0a <+11530>:	addps  %xmm7,%xmm11
   0x0000000000002d1d <+11549>:	addps  %xmm8,%xmm9
   0x0000000000002d5c <+11612>:	movaps 0x9c0(%r8),%xmm10
   0x0000000000002d64 <+11620>:	mulps  0x100(%rcx),%xmm10
   0x0000000000002da8 <+11688>:	addps  %xmm15,%xmm10

941
942	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
943	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
944	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002caa <+11434>:	movaps 0x910(%r8),%xmm9
   0x0000000000002cb2 <+11442>:	mulps  0x110(%rcx),%xmm9
   0x0000000000002cce <+11470>:	addps  %xmm11,%xmm9
   0x0000000000002d0e <+11534>:	movaps 0x950(%r8),%xmm7
   0x0000000000002d16 <+11542>:	mulps  0x110(%rcx),%xmm7
   0x0000000000002d21 <+11553>:	movaps 0x990(%r8),%xmm8
   0x0000000000002d29 <+11561>:	mulps  0x110(%rcx),%xmm8
   0x0000000000002d31 <+11569>:	addps  %xmm11,%xmm7
   0x0000000000002d6c <+11628>:	addps  %xmm9,%xmm8
   0x0000000000002dac <+11692>:	movaps 0x9d0(%r8),%xmm15
   0x0000000000002db4 <+11700>:	mulps  0x110(%rcx),%xmm15
   0x0000000000002dd0 <+11728>:	addps  %xmm10,%xmm15

945
946	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
947	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
948	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002cbe <+11454>:	movaps 0x920(%r8),%xmm10
   0x0000000000002cc6 <+11462>:	mulps  0x120(%rcx),%xmm10
   0x0000000000002ce2 <+11490>:	addps  %xmm9,%xmm10
   0x0000000000002cfa <+11514>:	movaps 0x960(%r8),%xmm10
   0x0000000000002d02 <+11522>:	mulps  0x120(%rcx),%xmm10
   0x0000000000002d35 <+11573>:	movaps 0x9a0(%r8),%xmm11
   0x0000000000002d3d <+11581>:	mulps  0x120(%rcx),%xmm11
   0x0000000000002d45 <+11589>:	addps  %xmm7,%xmm10
   0x0000000000002d80 <+11648>:	addps  %xmm8,%xmm11
   0x0000000000002dd4 <+11732>:	movaps 0x9e0(%r8),%xmm10
   0x0000000000002ddc <+11740>:	mulps  0x120(%rcx),%xmm10
   0x0000000000002df8 <+11768>:	addps  %xmm15,%xmm10

949
950	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
951	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
952	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002c27 <+11303>:	movaps 0x930(%r8),%xmm14
   0x0000000000002c2f <+11311>:	mulps  0x130(%rcx),%xmm14
   0x0000000000002c71 <+11377>:	movaps 0x9b0(%r8),%xmm12
   0x0000000000002c79 <+11385>:	mulps  0x130(%rcx),%xmm12
   0x0000000000002c85 <+11397>:	movaps 0x970(%r8),%xmm13
   0x0000000000002c8d <+11405>:	mulps  0x130(%rcx),%xmm13
   0x0000000000002cf6 <+11510>:	addps  %xmm10,%xmm14
   0x0000000000002d58 <+11608>:	addps  %xmm10,%xmm13
   0x0000000000002d94 <+11668>:	addps  %xmm11,%xmm12
   0x0000000000002d98 <+11672>:	movaps 0x9f0(%r8),%xmm11
   0x0000000000002da0 <+11680>:	mulps  0x130(%rcx),%xmm11
   0x0000000000002e0c <+11788>:	addps  %xmm10,%xmm11

953	      }
954
955	      /* A(3,3)*B(3,1) = C(3,1). */
956	      for (i = 0; i < 4; i++)
957	      {
958	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
959	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
960	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002d49 <+11593>:	movaps 0xa40(%r8),%xmm7
   0x0000000000002d51 <+11601>:	mulps  0x200(%rcx),%xmm7
   0x0000000000002d70 <+11632>:	movaps 0xa80(%r8),%xmm9
   0x0000000000002d78 <+11640>:	mulps  0x200(%rcx),%xmm9
   0x0000000000002d84 <+11652>:	movaps 0xa00(%r8),%xmm8
   0x0000000000002d8c <+11660>:	mulps  0x200(%rcx),%xmm8
   0x0000000000002dbc <+11708>:	addps  %xmm14,%xmm8
   0x0000000000002e48 <+11848>:	addps  %xmm13,%xmm7
   0x0000000000002e5c <+11868>:	addps  %xmm12,%xmm9
   0x0000000000002e9b <+11931>:	movaps 0xac0(%r8),%xmm13
   0x0000000000002ea3 <+11939>:	mulps  0x200(%rcx),%xmm13
   0x0000000000002ee6 <+12006>:	addps  %xmm11,%xmm13

961
962	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
963	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
964	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002dc0 <+11712>:	movaps 0xa10(%r8),%xmm14
   0x0000000000002dc8 <+11720>:	mulps  0x210(%rcx),%xmm14
   0x0000000000002de4 <+11748>:	addps  %xmm8,%xmm14
   0x0000000000002e4c <+11852>:	movaps 0xa50(%r8),%xmm13
   0x0000000000002e54 <+11860>:	mulps  0x210(%rcx),%xmm13
   0x0000000000002e60 <+11872>:	movaps 0xa90(%r8),%xmm12
   0x0000000000002e68 <+11880>:	mulps  0x210(%rcx),%xmm12
   0x0000000000002e70 <+11888>:	addps  %xmm7,%xmm13
   0x0000000000002e83 <+11907>:	addps  %xmm9,%xmm12
   0x0000000000002eea <+12010>:	movaps 0xad0(%r8),%xmm11
   0x0000000000002ef2 <+12018>:	mulps  0x210(%rcx),%xmm11
   0x0000000000002f0e <+12046>:	addps  %xmm13,%xmm11

965
966	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
967	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
968	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002de8 <+11752>:	movaps 0xa20(%r8),%xmm8
   0x0000000000002df0 <+11760>:	mulps  0x220(%rcx),%xmm8
   0x0000000000002e20 <+11808>:	addps  %xmm14,%xmm8
   0x0000000000002e24 <+11812>:	movaps 0xae0(%r8),%xmm14
   0x0000000000002e2c <+11820>:	mulps  0x220(%rcx),%xmm14
   0x0000000000002e74 <+11892>:	movaps 0xa60(%r8),%xmm7
   0x0000000000002e7c <+11900>:	mulps  0x220(%rcx),%xmm7
   0x0000000000002e87 <+11911>:	movaps 0xaa0(%r8),%xmm9
   0x0000000000002e8f <+11919>:	mulps  0x220(%rcx),%xmm9
   0x0000000000002e97 <+11927>:	addps  %xmm13,%xmm7
   0x0000000000002ebe <+11966>:	addps  %xmm12,%xmm9
   0x0000000000002f22 <+12066>:	addps  %xmm11,%xmm14

969
970	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
971	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
972	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002e10 <+11792>:	movaps 0xa30(%r8),%xmm10
   0x0000000000002e18 <+11800>:	mulps  0x230(%rcx),%xmm10
   0x0000000000002e34 <+11828>:	addps  %xmm8,%xmm10
   0x0000000000002e38 <+11832>:	movaps 0xa70(%r8),%xmm8
   0x0000000000002e40 <+11840>:	mulps  0x230(%rcx),%xmm8
   0x0000000000002eab <+11947>:	addps  %xmm7,%xmm8
   0x0000000000002eaf <+11951>:	movaps 0xab0(%r8),%xmm7
   0x0000000000002eb7 <+11959>:	mulps  0x230(%rcx),%xmm7
   0x0000000000002ed2 <+11986>:	addps  %xmm9,%xmm7
   0x0000000000002ed6 <+11990>:	movaps 0xaf0(%r8),%xmm9
   0x0000000000002ede <+11998>:	mulps  0x230(%rcx),%xmm9
   0x0000000000002f36 <+12086>:	addps  %xmm14,%xmm9

973	      }
974
975	      /* A(3,4)*B(4,1) = C(3,1). */
976	      for (i = 0; i < 4; i++)
977	      {
978	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
979	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
980	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002dfc <+11772>:	movaps 0xb00(%r8),%xmm15
   0x0000000000002e04 <+11780>:	mulps  0x300(%rcx),%xmm15
   0x0000000000002ec2 <+11970>:	movaps 0xb40(%r8),%xmm12
   0x0000000000002eca <+11978>:	mulps  0x300(%rcx),%xmm12
   0x0000000000002efa <+12026>:	addps  %xmm10,%xmm15
   0x0000000000002f3a <+12090>:	movaps 0xb80(%r8),%xmm14
   0x0000000000002f42 <+12098>:	mulps  0x300(%rcx),%xmm14
   0x0000000000002f62 <+12130>:	addps  %xmm8,%xmm12
   0x0000000000002f8e <+12174>:	movaps 0xbc0(%r8),%xmm12
   0x0000000000002f96 <+12182>:	mulps  0x300(%rcx),%xmm12
   0x0000000000002fb2 <+12210>:	addps  %xmm7,%xmm14
   0x0000000000002fc9 <+12233>:	addps  %xmm9,%xmm12

981
982	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
983	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
984	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002efe <+12030>:	movaps 0xb10(%r8),%xmm10
   0x0000000000002f06 <+12038>:	mulps  0x310(%rcx),%xmm10
   0x0000000000002f4a <+12106>:	addps  %xmm15,%xmm10
   0x0000000000002f66 <+12134>:	movaps 0xb50(%r8),%xmm8
   0x0000000000002f6e <+12142>:	mulps  0x310(%rcx),%xmm8
   0x0000000000002f8a <+12170>:	addps  %xmm12,%xmm8
   0x0000000000002fb6 <+12214>:	movaps 0xb90(%r8),%xmm7
   0x0000000000002fbe <+12222>:	mulps  0x310(%rcx),%xmm7
   0x0000000000002fcd <+12237>:	movaps 0xbd0(%r8),%xmm9
   0x0000000000002fd5 <+12245>:	mulps  0x310(%rcx),%xmm9
   0x0000000000002fdd <+12253>:	addps  %xmm14,%xmm7
   0x0000000000002ff4 <+12276>:	addps  %xmm12,%xmm9

985
986	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
987	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
988	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002f12 <+12050>:	movaps 0xb60(%r8),%xmm13
   0x0000000000002f1a <+12058>:	mulps  0x320(%rcx),%xmm13
   0x0000000000002f26 <+12070>:	movaps 0xb20(%r8),%xmm11
   0x0000000000002f2e <+12078>:	mulps  0x320(%rcx),%xmm11
   0x0000000000002f4e <+12110>:	addps  %xmm10,%xmm11
   0x0000000000002f7a <+12154>:	movaps 0xba0(%r8),%xmm11
   0x0000000000002f82 <+12162>:	mulps  0x320(%rcx),%xmm11
   0x0000000000002f9e <+12190>:	addps  %xmm8,%xmm13
   0x0000000000002fe1 <+12257>:	addps  %xmm7,%xmm11
   0x000000000000300c <+12300>:	movaps 0xbe0(%r8),%xmm11
   0x0000000000003014 <+12308>:	mulps  0x320(%rcx),%xmm11
   0x000000000000301c <+12316>:	addps  %xmm9,%xmm11

989
990	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
991	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
992	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000002f52 <+12114>:	movaps 0xb30(%r8),%xmm10
   0x0000000000002f5a <+12122>:	mulps  0x330(%rcx),%xmm10
   0x0000000000002f76 <+12150>:	addps  %xmm11,%xmm10
   0x0000000000002fa2 <+12194>:	movaps 0xb70(%r8),%xmm8
   0x0000000000002faa <+12202>:	mulps  0x330(%rcx),%xmm8
   0x0000000000002fc5 <+12229>:	addps  %xmm13,%xmm8
   0x0000000000002fe5 <+12261>:	movaps 0xbb0(%r8),%xmm7
   0x0000000000002fed <+12269>:	mulps  0x330(%rcx),%xmm7
   0x0000000000002ff8 <+12280>:	movaps 0xbf0(%r8),%xmm12
   0x0000000000003000 <+12288>:	mulps  0x330(%rcx),%xmm12
   0x0000000000003008 <+12296>:	addps  %xmm11,%xmm7
   0x0000000000003020 <+12320>:	addps  %xmm11,%xmm12

993	      }
994	    }
995
996	    /* Store C(3,1) block. */
997	    for (i = 0; i < 4; i++)
998	    {
999	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000003024 <+12324>:	mulps  %xmm6,%xmm10
   0x0000000000003040 <+12352>:	mulps  %xmm6,%xmm8
   0x0000000000003058 <+12376>:	mulps  %xmm6,%xmm7
   0x000000000000306e <+12398>:	mulps  %xmm6,%xmm12

1000	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_31]), C_row[i]);
   0x0000000000003028 <+12328>:	addps  0x200(%r9),%xmm10
   0x0000000000003044 <+12356>:	addps  0x210(%r9),%xmm8
   0x000000000000305b <+12379>:	addps  0x220(%r9),%xmm7
   0x0000000000003072 <+12402>:	addps  0x230(%r9),%xmm12

1001	      _mm_store_ps(&C[i*4+C_OFFSET_31], C_row[i]);
   0x0000000000003034 <+12340>:	movaps %xmm10,0x200(%r9)
   0x000000000000304c <+12364>:	movaps %xmm8,0x210(%r9)
   0x0000000000003063 <+12387>:	movaps %xmm7,0x220(%r9)
   0x000000000000307a <+12410>:	movaps %xmm12,0x230(%r9)

1002	    }
1003
1004	    /* Reset C(3,2) matrix accumulators */
1005	    C_row[0] = _mm_setzero_ps();
   0x000000000000303c <+12348>:	xorps  %xmm10,%xmm10

1006	    C_row[1] = _mm_setzero_ps();
   0x0000000000003054 <+12372>:	xorps  %xmm8,%xmm8

1007	    C_row[2] = _mm_setzero_ps();
   0x000000000000306b <+12395>:	xorps  %xmm7,%xmm7

1008	    C_row[3] = _mm_setzero_ps();
   0x0000000000003082 <+12418>:	xorps  %xmm12,%xmm12

1009
1010	    if (norm[8]*norm[17] >= tolerance &&
   0x0000000000003030 <+12336>:	movaps %xmm4,%xmm9
   0x0000000000003086 <+12422>:	mulss  %xmm0,%xmm9
   0x000000000000308b <+12427>:	comiss %xmm1,%xmm9
   0x000000000000308f <+12431>:	jb     0x3599 <stream_kernel+13721>

1011	        norm[9]*norm[21] >= tolerance &&
   0x0000000000003095 <+12437>:	movss  0x3c(%rdx,%rsi,1),%xmm9
   0x000000000000309c <+12444>:	mulss  0x6c(%rdx,%rsi,1),%xmm9
   0x00000000000030a3 <+12451>:	comiss %xmm1,%xmm9
   0x00000000000030a7 <+12455>:	jb     0x3599 <stream_kernel+13721>

1012	        norm[10]*norm[25] >= tolerance &&
   0x00000000000030ad <+12461>:	movss  0x40(%rdx,%rsi,1),%xmm9
   0x00000000000030b4 <+12468>:	mulss  0x7c(%rdx,%rsi,1),%xmm9
   0x00000000000030bb <+12475>:	comiss %xmm1,%xmm9
   0x00000000000030bf <+12479>:	jb     0x3599 <stream_kernel+13721>

1013	        norm[11]*norm[29] >= tolerance)
   0x00000000000030c5 <+12485>:	movss  0x44(%rdx,%rsi,1),%xmm9
   0x00000000000030cc <+12492>:	mulss  0x8c(%rdx,%rsi,1),%xmm9
   0x00000000000030d6 <+12502>:	comiss %xmm1,%xmm9
   0x00000000000030da <+12506>:	jb     0x3599 <stream_kernel+13721>

1014	    {
1015	      /* A(3,1)*B(1,2) = C(3,2). */
1016	      for (i = 0; i < 4; i++)
1017	      {
1018	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
1019	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
1020	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000030e0 <+12512>:	movaps 0x800(%r8),%xmm10
   0x0000000000003100 <+12544>:	movaps 0x840(%r8),%xmm13
   0x0000000000003118 <+12568>:	movaps 0x880(%r8),%xmm8
   0x0000000000003120 <+12576>:	mulps  0x40(%rcx),%xmm10
   0x0000000000003134 <+12596>:	mulps  0x40(%rcx),%xmm13
   0x0000000000003142 <+12610>:	mulps  0x40(%rcx),%xmm8
   0x000000000000318b <+12683>:	movaps 0x8c0(%r8),%xmm13
   0x0000000000003193 <+12691>:	mulps  0x40(%rcx),%xmm13

1021
1022	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
1023	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
1024	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000030e8 <+12520>:	movaps 0x810(%r8),%xmm11
   0x0000000000003108 <+12552>:	movaps 0x850(%r8),%xmm14
   0x0000000000003125 <+12581>:	mulps  0x50(%rcx),%xmm11
   0x0000000000003139 <+12601>:	mulps  0x50(%rcx),%xmm14
   0x0000000000003147 <+12615>:	movaps 0x8d0(%r8),%xmm15
   0x000000000000314f <+12623>:	addps  %xmm10,%xmm11
   0x0000000000003153 <+12627>:	mulps  0x50(%rcx),%xmm15
   0x0000000000003164 <+12644>:	movaps 0x890(%r8),%xmm11
   0x0000000000003171 <+12657>:	mulps  0x50(%rcx),%xmm11
   0x0000000000003187 <+12679>:	addps  %xmm13,%xmm14
   0x00000000000031bd <+12733>:	addps  %xmm8,%xmm11
   0x00000000000031f6 <+12790>:	addps  %xmm13,%xmm15

1025
1026	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
1027	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
1028	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000030f0 <+12528>:	movaps 0x820(%r8),%xmm12
   0x000000000000312a <+12586>:	mulps  0x60(%rcx),%xmm12
   0x0000000000003158 <+12632>:	movaps 0x860(%r8),%xmm10
   0x0000000000003160 <+12640>:	addps  %xmm11,%xmm12
   0x000000000000316c <+12652>:	mulps  0x60(%rcx),%xmm10
   0x000000000000317a <+12666>:	movaps 0x8a0(%r8),%xmm12
   0x0000000000003182 <+12674>:	mulps  0x60(%rcx),%xmm12
   0x0000000000003198 <+12696>:	addps  %xmm14,%xmm10
   0x00000000000031b0 <+12720>:	movaps 0x8e0(%r8),%xmm10
   0x00000000000031b8 <+12728>:	mulps  0x60(%rcx),%xmm10
   0x00000000000031ce <+12750>:	addps  %xmm11,%xmm12
   0x000000000000320a <+12810>:	addps  %xmm15,%xmm10

1029
1030	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
1031	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
1032	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000030f8 <+12536>:	movaps 0x830(%r8),%xmm9
   0x0000000000003110 <+12560>:	movaps 0x870(%r8),%xmm7
   0x000000000000312f <+12591>:	mulps  0x70(%rcx),%xmm9
   0x000000000000313e <+12606>:	mulps  0x70(%rcx),%xmm7
   0x0000000000003176 <+12662>:	addps  %xmm12,%xmm9
   0x00000000000031ac <+12716>:	addps  %xmm10,%xmm7
   0x00000000000031c1 <+12737>:	movaps 0x8b0(%r8),%xmm8
   0x00000000000031c9 <+12745>:	mulps  0x70(%rcx),%xmm8
   0x00000000000031e2 <+12770>:	addps  %xmm12,%xmm8
   0x000000000000320e <+12814>:	movaps 0x8f0(%r8),%xmm15
   0x0000000000003216 <+12822>:	mulps  0x70(%rcx),%xmm15
   0x000000000000322f <+12847>:	addps  %xmm10,%xmm15

1033	      }
1034
1035	      /* A(3,2)*B(2,2) = C(3,2). */
1036	      for (i = 0; i < 4; i++)
1037	      {
1038	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1039	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
1040	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000031d2 <+12754>:	movaps 0x900(%r8),%xmm11
   0x00000000000031da <+12762>:	mulps  0x140(%rcx),%xmm11
   0x000000000000321b <+12827>:	addps  %xmm9,%xmm11
   0x0000000000003247 <+12871>:	movaps 0x940(%r8),%xmm11
   0x000000000000324f <+12879>:	mulps  0x140(%rcx),%xmm11
   0x000000000000325b <+12891>:	movaps 0x980(%r8),%xmm9
   0x0000000000003263 <+12899>:	mulps  0x140(%rcx),%xmm9
   0x000000000000327f <+12927>:	addps  %xmm7,%xmm11
   0x0000000000003292 <+12946>:	addps  %xmm8,%xmm9
   0x00000000000032d1 <+13009>:	movaps 0x9c0(%r8),%xmm10
   0x00000000000032d9 <+13017>:	mulps  0x140(%rcx),%xmm10
   0x000000000000331d <+13085>:	addps  %xmm15,%xmm10

1041
1042	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1043	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
1044	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000321f <+12831>:	movaps 0x910(%r8),%xmm9
   0x0000000000003227 <+12839>:	mulps  0x150(%rcx),%xmm9
   0x0000000000003243 <+12867>:	addps  %xmm11,%xmm9
   0x0000000000003283 <+12931>:	movaps 0x950(%r8),%xmm7
   0x000000000000328b <+12939>:	mulps  0x150(%rcx),%xmm7
   0x0000000000003296 <+12950>:	movaps 0x990(%r8),%xmm8
   0x000000000000329e <+12958>:	mulps  0x150(%rcx),%xmm8
   0x00000000000032a6 <+12966>:	addps  %xmm11,%xmm7
   0x00000000000032e1 <+13025>:	addps  %xmm9,%xmm8
   0x0000000000003321 <+13089>:	movaps 0x9d0(%r8),%xmm15
   0x0000000000003329 <+13097>:	mulps  0x150(%rcx),%xmm15
   0x0000000000003345 <+13125>:	addps  %xmm10,%xmm15

1045
1046	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1047	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
1048	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003233 <+12851>:	movaps 0x920(%r8),%xmm10
   0x000000000000323b <+12859>:	mulps  0x160(%rcx),%xmm10
   0x0000000000003257 <+12887>:	addps  %xmm9,%xmm10
   0x000000000000326f <+12911>:	movaps 0x960(%r8),%xmm10
   0x0000000000003277 <+12919>:	mulps  0x160(%rcx),%xmm10
   0x00000000000032aa <+12970>:	movaps 0x9a0(%r8),%xmm11
   0x00000000000032b2 <+12978>:	mulps  0x160(%rcx),%xmm11
   0x00000000000032ba <+12986>:	addps  %xmm7,%xmm10
   0x00000000000032f5 <+13045>:	addps  %xmm8,%xmm11
   0x0000000000003349 <+13129>:	movaps 0x9e0(%r8),%xmm10
   0x0000000000003351 <+13137>:	mulps  0x160(%rcx),%xmm10
   0x000000000000336d <+13165>:	addps  %xmm15,%xmm10

1049
1050	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1051	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
1052	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000319c <+12700>:	movaps 0x930(%r8),%xmm14
   0x00000000000031a4 <+12708>:	mulps  0x170(%rcx),%xmm14
   0x00000000000031e6 <+12774>:	movaps 0x9b0(%r8),%xmm12
   0x00000000000031ee <+12782>:	mulps  0x170(%rcx),%xmm12
   0x00000000000031fa <+12794>:	movaps 0x970(%r8),%xmm13
   0x0000000000003202 <+12802>:	mulps  0x170(%rcx),%xmm13
   0x000000000000326b <+12907>:	addps  %xmm10,%xmm14
   0x00000000000032cd <+13005>:	addps  %xmm10,%xmm13
   0x0000000000003309 <+13065>:	addps  %xmm11,%xmm12
   0x000000000000330d <+13069>:	movaps 0x9f0(%r8),%xmm11
   0x0000000000003315 <+13077>:	mulps  0x170(%rcx),%xmm11
   0x0000000000003381 <+13185>:	addps  %xmm10,%xmm11

1053	      }
1054
1055	      /* A(3,3)*B(3,2) = C(3,2). */
1056	      for (i = 0; i < 4; i++)
1057	      {
1058	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1059	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
1060	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000032be <+12990>:	movaps 0xa40(%r8),%xmm7
   0x00000000000032c6 <+12998>:	mulps  0x240(%rcx),%xmm7
   0x00000000000032e5 <+13029>:	movaps 0xa80(%r8),%xmm9
   0x00000000000032ed <+13037>:	mulps  0x240(%rcx),%xmm9
   0x00000000000032f9 <+13049>:	movaps 0xa00(%r8),%xmm8
   0x0000000000003301 <+13057>:	mulps  0x240(%rcx),%xmm8
   0x0000000000003331 <+13105>:	addps  %xmm14,%xmm8
   0x00000000000033bd <+13245>:	addps  %xmm13,%xmm7
   0x00000000000033d1 <+13265>:	addps  %xmm12,%xmm9
   0x0000000000003410 <+13328>:	movaps 0xac0(%r8),%xmm13
   0x0000000000003418 <+13336>:	mulps  0x240(%rcx),%xmm13
   0x000000000000345b <+13403>:	addps  %xmm11,%xmm13

1061
1062	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1063	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
1064	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003335 <+13109>:	movaps 0xa10(%r8),%xmm14
   0x000000000000333d <+13117>:	mulps  0x250(%rcx),%xmm14
   0x0000000000003359 <+13145>:	addps  %xmm8,%xmm14
   0x00000000000033c1 <+13249>:	movaps 0xa50(%r8),%xmm13
   0x00000000000033c9 <+13257>:	mulps  0x250(%rcx),%xmm13
   0x00000000000033d5 <+13269>:	movaps 0xa90(%r8),%xmm12
   0x00000000000033dd <+13277>:	mulps  0x250(%rcx),%xmm12
   0x00000000000033e5 <+13285>:	addps  %xmm7,%xmm13
   0x00000000000033f8 <+13304>:	addps  %xmm9,%xmm12
   0x000000000000345f <+13407>:	movaps 0xad0(%r8),%xmm11
   0x0000000000003467 <+13415>:	mulps  0x250(%rcx),%xmm11
   0x0000000000003483 <+13443>:	addps  %xmm13,%xmm11

1065
1066	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1067	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
1068	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000335d <+13149>:	movaps 0xa20(%r8),%xmm8
   0x0000000000003365 <+13157>:	mulps  0x260(%rcx),%xmm8
   0x0000000000003395 <+13205>:	addps  %xmm14,%xmm8
   0x0000000000003399 <+13209>:	movaps 0xae0(%r8),%xmm14
   0x00000000000033a1 <+13217>:	mulps  0x260(%rcx),%xmm14
   0x00000000000033e9 <+13289>:	movaps 0xa60(%r8),%xmm7
   0x00000000000033f1 <+13297>:	mulps  0x260(%rcx),%xmm7
   0x00000000000033fc <+13308>:	movaps 0xaa0(%r8),%xmm9
   0x0000000000003404 <+13316>:	mulps  0x260(%rcx),%xmm9
   0x000000000000340c <+13324>:	addps  %xmm13,%xmm7
   0x0000000000003433 <+13363>:	addps  %xmm12,%xmm9
   0x0000000000003497 <+13463>:	addps  %xmm11,%xmm14

1069
1070	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1071	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
1072	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003385 <+13189>:	movaps 0xa30(%r8),%xmm10
   0x000000000000338d <+13197>:	mulps  0x270(%rcx),%xmm10
   0x00000000000033a9 <+13225>:	addps  %xmm8,%xmm10
   0x00000000000033ad <+13229>:	movaps 0xa70(%r8),%xmm8
   0x00000000000033b5 <+13237>:	mulps  0x270(%rcx),%xmm8
   0x0000000000003420 <+13344>:	addps  %xmm7,%xmm8
   0x0000000000003424 <+13348>:	movaps 0xab0(%r8),%xmm7
   0x000000000000342c <+13356>:	mulps  0x270(%rcx),%xmm7
   0x0000000000003447 <+13383>:	addps  %xmm9,%xmm7
   0x000000000000344b <+13387>:	movaps 0xaf0(%r8),%xmm9
   0x0000000000003453 <+13395>:	mulps  0x270(%rcx),%xmm9
   0x00000000000034ab <+13483>:	addps  %xmm14,%xmm9

1073	      }
1074
1075	      /* A(3,4)*B(4,2) = C(3,2). */
1076	      for (i = 0; i < 4; i++)
1077	      {
1078	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1079	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
1080	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003371 <+13169>:	movaps 0xb00(%r8),%xmm15
   0x0000000000003379 <+13177>:	mulps  0x340(%rcx),%xmm15
   0x0000000000003437 <+13367>:	movaps 0xb40(%r8),%xmm12
   0x000000000000343f <+13375>:	mulps  0x340(%rcx),%xmm12
   0x000000000000346f <+13423>:	addps  %xmm10,%xmm15
   0x00000000000034af <+13487>:	movaps 0xb80(%r8),%xmm14
   0x00000000000034b7 <+13495>:	mulps  0x340(%rcx),%xmm14
   0x00000000000034d7 <+13527>:	addps  %xmm8,%xmm12
   0x0000000000003503 <+13571>:	movaps 0xbc0(%r8),%xmm12
   0x000000000000350b <+13579>:	mulps  0x340(%rcx),%xmm12
   0x0000000000003527 <+13607>:	addps  %xmm7,%xmm14
   0x000000000000353e <+13630>:	addps  %xmm9,%xmm12

1081
1082	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1083	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
1084	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003473 <+13427>:	movaps 0xb10(%r8),%xmm10
   0x000000000000347b <+13435>:	mulps  0x350(%rcx),%xmm10
   0x00000000000034bf <+13503>:	addps  %xmm15,%xmm10
   0x00000000000034db <+13531>:	movaps 0xb50(%r8),%xmm8
   0x00000000000034e3 <+13539>:	mulps  0x350(%rcx),%xmm8
   0x00000000000034ff <+13567>:	addps  %xmm12,%xmm8
   0x000000000000352b <+13611>:	movaps 0xb90(%r8),%xmm7
   0x0000000000003533 <+13619>:	mulps  0x350(%rcx),%xmm7
   0x0000000000003542 <+13634>:	movaps 0xbd0(%r8),%xmm9
   0x000000000000354a <+13642>:	mulps  0x350(%rcx),%xmm9
   0x0000000000003552 <+13650>:	addps  %xmm14,%xmm7
   0x0000000000003569 <+13673>:	addps  %xmm12,%xmm9

1085
1086	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1087	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
1088	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003487 <+13447>:	movaps 0xb60(%r8),%xmm13
   0x000000000000348f <+13455>:	mulps  0x360(%rcx),%xmm13
   0x000000000000349b <+13467>:	movaps 0xb20(%r8),%xmm11
   0x00000000000034a3 <+13475>:	mulps  0x360(%rcx),%xmm11
   0x00000000000034c3 <+13507>:	addps  %xmm10,%xmm11
   0x00000000000034ef <+13551>:	movaps 0xba0(%r8),%xmm11
   0x00000000000034f7 <+13559>:	mulps  0x360(%rcx),%xmm11
   0x0000000000003513 <+13587>:	addps  %xmm8,%xmm13
   0x0000000000003556 <+13654>:	addps  %xmm7,%xmm11
   0x0000000000003581 <+13697>:	movaps 0xbe0(%r8),%xmm11
   0x0000000000003589 <+13705>:	mulps  0x360(%rcx),%xmm11
   0x0000000000003591 <+13713>:	addps  %xmm9,%xmm11

1089
1090	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1091	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
1092	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000034c7 <+13511>:	movaps 0xb30(%r8),%xmm10
   0x00000000000034cf <+13519>:	mulps  0x370(%rcx),%xmm10
   0x00000000000034eb <+13547>:	addps  %xmm11,%xmm10
   0x0000000000003517 <+13591>:	movaps 0xb70(%r8),%xmm8
   0x000000000000351f <+13599>:	mulps  0x370(%rcx),%xmm8
   0x000000000000353a <+13626>:	addps  %xmm13,%xmm8
   0x000000000000355a <+13658>:	movaps 0xbb0(%r8),%xmm7
   0x0000000000003562 <+13666>:	mulps  0x370(%rcx),%xmm7
   0x000000000000356d <+13677>:	movaps 0xbf0(%r8),%xmm12
   0x0000000000003575 <+13685>:	mulps  0x370(%rcx),%xmm12
   0x000000000000357d <+13693>:	addps  %xmm11,%xmm7
   0x0000000000003595 <+13717>:	addps  %xmm11,%xmm12

1093	      }
1094	    }
1095
1096	    /* Store C(3,2) block. */
1097	    for (i = 0; i < 4; i++)
1098	    {
1099	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000003599 <+13721>:	mulps  %xmm6,%xmm10
   0x00000000000035b5 <+13749>:	mulps  %xmm6,%xmm8
   0x00000000000035cd <+13773>:	mulps  %xmm6,%xmm7
   0x00000000000035e3 <+13795>:	mulps  %xmm6,%xmm12

1100	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_32]), C_row[i]);
   0x000000000000359d <+13725>:	addps  0x240(%r9),%xmm10
   0x00000000000035b9 <+13753>:	addps  0x250(%r9),%xmm8
   0x00000000000035d0 <+13776>:	addps  0x260(%r9),%xmm7
   0x00000000000035e7 <+13799>:	addps  0x270(%r9),%xmm12

1101	      _mm_store_ps(&C[i*4+C_OFFSET_32], C_row[i]);
   0x00000000000035a9 <+13737>:	movaps %xmm10,0x240(%r9)
   0x00000000000035c1 <+13761>:	movaps %xmm8,0x250(%r9)
   0x00000000000035d8 <+13784>:	movaps %xmm7,0x260(%r9)
   0x00000000000035ef <+13807>:	movaps %xmm12,0x270(%r9)

1102	    }
1103
1104	    /* Reset C(3,3) matrix accumulators */
1105	    C_row[0] = _mm_setzero_ps();
   0x00000000000035b1 <+13745>:	xorps  %xmm10,%xmm10

1106	    C_row[1] = _mm_setzero_ps();
   0x00000000000035c9 <+13769>:	xorps  %xmm8,%xmm8

1107	    C_row[2] = _mm_setzero_ps();
   0x00000000000035e0 <+13792>:	xorps  %xmm7,%xmm7

1108	    C_row[3] = _mm_setzero_ps();
   0x00000000000035f7 <+13815>:	xorps  %xmm12,%xmm12

1109
1110	    if (norm[8]*norm[18] >= tolerance &&
   0x00000000000035a5 <+13733>:	movaps %xmm3,%xmm9
   0x00000000000035fb <+13819>:	mulss  %xmm0,%xmm9
   0x0000000000003600 <+13824>:	comiss %xmm1,%xmm9
   0x0000000000003604 <+13828>:	jb     0x3b41 <stream_kernel+15169>

1111	        norm[9]*norm[22] >= tolerance &&
   0x000000000000360a <+13834>:	movss  0x3c(%rdx,%rsi,1),%xmm9
   0x0000000000003611 <+13841>:	mulss  0x70(%rdx,%rsi,1),%xmm9
   0x0000000000003618 <+13848>:	comiss %xmm1,%xmm9
   0x000000000000361c <+13852>:	jb     0x3b41 <stream_kernel+15169>

1112	        norm[10]*norm[26] >= tolerance &&
   0x0000000000003622 <+13858>:	movss  0x40(%rdx,%rsi,1),%xmm9
   0x0000000000003629 <+13865>:	mulss  0x80(%rdx,%rsi,1),%xmm9
   0x0000000000003633 <+13875>:	comiss %xmm1,%xmm9
   0x0000000000003637 <+13879>:	jb     0x3b41 <stream_kernel+15169>

1113	        norm[11]*norm[30] >= tolerance)
   0x000000000000363d <+13885>:	movss  0x44(%rdx,%rsi,1),%xmm9
   0x0000000000003644 <+13892>:	mulss  0x90(%rdx,%rsi,1),%xmm9
   0x000000000000364e <+13902>:	comiss %xmm1,%xmm9
   0x0000000000003652 <+13906>:	jb     0x3b41 <stream_kernel+15169>

1114	    {
1115	      /* A(3,1)*B(1,3) = C(3,3). */
1116	      for (i = 0; i < 4; i++)
1117	      {
1118	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
1119	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
1120	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003658 <+13912>:	movaps 0x800(%r8),%xmm10
   0x0000000000003678 <+13944>:	movaps 0x840(%r8),%xmm13
   0x0000000000003690 <+13968>:	movaps 0x880(%r8),%xmm8
   0x0000000000003698 <+13976>:	mulps  0x80(%rcx),%xmm10
   0x00000000000036b8 <+14008>:	mulps  0x80(%rcx),%xmm13
   0x00000000000036cf <+14031>:	mulps  0x80(%rcx),%xmm8
   0x0000000000003727 <+14119>:	movaps 0x8c0(%r8),%xmm13
   0x000000000000372f <+14127>:	mulps  0x80(%rcx),%xmm13

1121
1122	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
1123	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
1124	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003660 <+13920>:	movaps 0x810(%r8),%xmm11
   0x0000000000003680 <+13952>:	movaps 0x850(%r8),%xmm14
   0x00000000000036a0 <+13984>:	mulps  0x90(%rcx),%xmm11
   0x00000000000036c0 <+14016>:	mulps  0x90(%rcx),%xmm14
   0x00000000000036d7 <+14039>:	movaps 0x8d0(%r8),%xmm15
   0x00000000000036df <+14047>:	addps  %xmm10,%xmm11
   0x00000000000036e3 <+14051>:	mulps  0x90(%rcx),%xmm15
   0x00000000000036f7 <+14071>:	movaps 0x890(%r8),%xmm11
   0x0000000000003707 <+14087>:	mulps  0x90(%rcx),%xmm11
   0x0000000000003723 <+14115>:	addps  %xmm13,%xmm14
   0x000000000000375f <+14175>:	addps  %xmm8,%xmm11
   0x000000000000379b <+14235>:	addps  %xmm13,%xmm15

1125
1126	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
1127	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
1128	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003668 <+13928>:	movaps 0x820(%r8),%xmm12
   0x00000000000036a8 <+13992>:	mulps  0xa0(%rcx),%xmm12
   0x00000000000036eb <+14059>:	movaps 0x860(%r8),%xmm10
   0x00000000000036f3 <+14067>:	addps  %xmm11,%xmm12
   0x00000000000036ff <+14079>:	mulps  0xa0(%rcx),%xmm10
   0x0000000000003713 <+14099>:	movaps 0x8a0(%r8),%xmm12
   0x000000000000371b <+14107>:	mulps  0xa0(%rcx),%xmm12
   0x0000000000003737 <+14135>:	addps  %xmm14,%xmm10
   0x000000000000374f <+14159>:	movaps 0x8e0(%r8),%xmm10
   0x0000000000003757 <+14167>:	mulps  0xa0(%rcx),%xmm10
   0x0000000000003773 <+14195>:	addps  %xmm11,%xmm12
   0x00000000000037af <+14255>:	addps  %xmm15,%xmm10

1129
1130	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
1131	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
1132	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003670 <+13936>:	movaps 0x830(%r8),%xmm9
   0x0000000000003688 <+13960>:	movaps 0x870(%r8),%xmm7
   0x00000000000036b0 <+14000>:	mulps  0xb0(%rcx),%xmm9
   0x00000000000036c8 <+14024>:	mulps  0xb0(%rcx),%xmm7
   0x000000000000370f <+14095>:	addps  %xmm12,%xmm9
   0x000000000000374b <+14155>:	addps  %xmm10,%xmm7
   0x0000000000003763 <+14179>:	movaps 0x8b0(%r8),%xmm8
   0x000000000000376b <+14187>:	mulps  0xb0(%rcx),%xmm8
   0x0000000000003787 <+14215>:	addps  %xmm12,%xmm8
   0x00000000000037b3 <+14259>:	movaps 0x8f0(%r8),%xmm15
   0x00000000000037bb <+14267>:	mulps  0xb0(%rcx),%xmm15
   0x00000000000037d7 <+14295>:	addps  %xmm10,%xmm15

1133	      }
1134
1135	      /* A(3,2)*B(2,3) = C(3,3). */
1136	      for (i = 0; i < 4; i++)
1137	      {
1138	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1139	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
1140	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003777 <+14199>:	movaps 0x900(%r8),%xmm11
   0x000000000000377f <+14207>:	mulps  0x180(%rcx),%xmm11
   0x00000000000037c3 <+14275>:	addps  %xmm9,%xmm11
   0x00000000000037ef <+14319>:	movaps 0x940(%r8),%xmm11
   0x00000000000037f7 <+14327>:	mulps  0x180(%rcx),%xmm11
   0x0000000000003803 <+14339>:	movaps 0x980(%r8),%xmm9
   0x000000000000380b <+14347>:	mulps  0x180(%rcx),%xmm9
   0x0000000000003827 <+14375>:	addps  %xmm7,%xmm11
   0x000000000000383a <+14394>:	addps  %xmm8,%xmm9
   0x0000000000003879 <+14457>:	movaps 0x9c0(%r8),%xmm10
   0x0000000000003881 <+14465>:	mulps  0x180(%rcx),%xmm10
   0x00000000000038c5 <+14533>:	addps  %xmm15,%xmm10

1141
1142	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1143	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
1144	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000037c7 <+14279>:	movaps 0x910(%r8),%xmm9
   0x00000000000037cf <+14287>:	mulps  0x190(%rcx),%xmm9
   0x00000000000037eb <+14315>:	addps  %xmm11,%xmm9
   0x000000000000382b <+14379>:	movaps 0x950(%r8),%xmm7
   0x0000000000003833 <+14387>:	mulps  0x190(%rcx),%xmm7
   0x000000000000383e <+14398>:	movaps 0x990(%r8),%xmm8
   0x0000000000003846 <+14406>:	mulps  0x190(%rcx),%xmm8
   0x000000000000384e <+14414>:	addps  %xmm11,%xmm7
   0x0000000000003889 <+14473>:	addps  %xmm9,%xmm8
   0x00000000000038c9 <+14537>:	movaps 0x9d0(%r8),%xmm15
   0x00000000000038d1 <+14545>:	mulps  0x190(%rcx),%xmm15
   0x00000000000038ed <+14573>:	addps  %xmm10,%xmm15

1145
1146	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1147	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
1148	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000037db <+14299>:	movaps 0x920(%r8),%xmm10
   0x00000000000037e3 <+14307>:	mulps  0x1a0(%rcx),%xmm10
   0x00000000000037ff <+14335>:	addps  %xmm9,%xmm10
   0x0000000000003817 <+14359>:	movaps 0x960(%r8),%xmm10
   0x000000000000381f <+14367>:	mulps  0x1a0(%rcx),%xmm10
   0x0000000000003852 <+14418>:	movaps 0x9a0(%r8),%xmm11
   0x000000000000385a <+14426>:	mulps  0x1a0(%rcx),%xmm11
   0x0000000000003862 <+14434>:	addps  %xmm7,%xmm10
   0x000000000000389d <+14493>:	addps  %xmm8,%xmm11
   0x00000000000038f1 <+14577>:	movaps 0x9e0(%r8),%xmm10
   0x00000000000038f9 <+14585>:	mulps  0x1a0(%rcx),%xmm10
   0x0000000000003915 <+14613>:	addps  %xmm15,%xmm10

1149
1150	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1151	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
1152	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000373b <+14139>:	movaps 0x930(%r8),%xmm14
   0x0000000000003743 <+14147>:	mulps  0x1b0(%rcx),%xmm14
   0x000000000000378b <+14219>:	movaps 0x9b0(%r8),%xmm12
   0x0000000000003793 <+14227>:	mulps  0x1b0(%rcx),%xmm12
   0x000000000000379f <+14239>:	movaps 0x970(%r8),%xmm13
   0x00000000000037a7 <+14247>:	mulps  0x1b0(%rcx),%xmm13
   0x0000000000003813 <+14355>:	addps  %xmm10,%xmm14
   0x0000000000003875 <+14453>:	addps  %xmm10,%xmm13
   0x00000000000038b1 <+14513>:	addps  %xmm11,%xmm12
   0x00000000000038b5 <+14517>:	movaps 0x9f0(%r8),%xmm11
   0x00000000000038bd <+14525>:	mulps  0x1b0(%rcx),%xmm11
   0x0000000000003929 <+14633>:	addps  %xmm10,%xmm11

1153	      }
1154
1155	      /* A(3,3)*B(3,3) = C(3,3). */
1156	      for (i = 0; i < 4; i++)
1157	      {
1158	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1159	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
1160	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003866 <+14438>:	movaps 0xa40(%r8),%xmm7
   0x000000000000386e <+14446>:	mulps  0x280(%rcx),%xmm7
   0x000000000000388d <+14477>:	movaps 0xa80(%r8),%xmm9
   0x0000000000003895 <+14485>:	mulps  0x280(%rcx),%xmm9
   0x00000000000038a1 <+14497>:	movaps 0xa00(%r8),%xmm8
   0x00000000000038a9 <+14505>:	mulps  0x280(%rcx),%xmm8
   0x00000000000038d9 <+14553>:	addps  %xmm14,%xmm8
   0x0000000000003965 <+14693>:	addps  %xmm13,%xmm7
   0x0000000000003979 <+14713>:	addps  %xmm12,%xmm9
   0x00000000000039b8 <+14776>:	movaps 0xac0(%r8),%xmm13
   0x00000000000039c0 <+14784>:	mulps  0x280(%rcx),%xmm13
   0x0000000000003a03 <+14851>:	addps  %xmm11,%xmm13

1161
1162	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1163	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
1164	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000038dd <+14557>:	movaps 0xa10(%r8),%xmm14
   0x00000000000038e5 <+14565>:	mulps  0x290(%rcx),%xmm14
   0x0000000000003901 <+14593>:	addps  %xmm8,%xmm14
   0x0000000000003969 <+14697>:	movaps 0xa50(%r8),%xmm13
   0x0000000000003971 <+14705>:	mulps  0x290(%rcx),%xmm13
   0x000000000000397d <+14717>:	movaps 0xa90(%r8),%xmm12
   0x0000000000003985 <+14725>:	mulps  0x290(%rcx),%xmm12
   0x000000000000398d <+14733>:	addps  %xmm7,%xmm13
   0x00000000000039a0 <+14752>:	addps  %xmm9,%xmm12
   0x0000000000003a07 <+14855>:	movaps 0xad0(%r8),%xmm11
   0x0000000000003a0f <+14863>:	mulps  0x290(%rcx),%xmm11
   0x0000000000003a2b <+14891>:	addps  %xmm13,%xmm11

1165
1166	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1167	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
1168	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003905 <+14597>:	movaps 0xa20(%r8),%xmm8
   0x000000000000390d <+14605>:	mulps  0x2a0(%rcx),%xmm8
   0x000000000000393d <+14653>:	addps  %xmm14,%xmm8
   0x0000000000003941 <+14657>:	movaps 0xae0(%r8),%xmm14
   0x0000000000003949 <+14665>:	mulps  0x2a0(%rcx),%xmm14
   0x0000000000003991 <+14737>:	movaps 0xa60(%r8),%xmm7
   0x0000000000003999 <+14745>:	mulps  0x2a0(%rcx),%xmm7
   0x00000000000039a4 <+14756>:	movaps 0xaa0(%r8),%xmm9
   0x00000000000039ac <+14764>:	mulps  0x2a0(%rcx),%xmm9
   0x00000000000039b4 <+14772>:	addps  %xmm13,%xmm7
   0x00000000000039db <+14811>:	addps  %xmm12,%xmm9
   0x0000000000003a3f <+14911>:	addps  %xmm11,%xmm14

1169
1170	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1171	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
1172	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000392d <+14637>:	movaps 0xa30(%r8),%xmm10
   0x0000000000003935 <+14645>:	mulps  0x2b0(%rcx),%xmm10
   0x0000000000003951 <+14673>:	addps  %xmm8,%xmm10
   0x0000000000003955 <+14677>:	movaps 0xa70(%r8),%xmm8
   0x000000000000395d <+14685>:	mulps  0x2b0(%rcx),%xmm8
   0x00000000000039c8 <+14792>:	addps  %xmm7,%xmm8
   0x00000000000039cc <+14796>:	movaps 0xab0(%r8),%xmm7
   0x00000000000039d4 <+14804>:	mulps  0x2b0(%rcx),%xmm7
   0x00000000000039ef <+14831>:	addps  %xmm9,%xmm7
   0x00000000000039f3 <+14835>:	movaps 0xaf0(%r8),%xmm9
   0x00000000000039fb <+14843>:	mulps  0x2b0(%rcx),%xmm9
   0x0000000000003a53 <+14931>:	addps  %xmm14,%xmm9

1173	      }
1174
1175	      /* A(3,4)*B(4,3) = C(3,3). */
1176	      for (i = 0; i < 4; i++)
1177	      {
1178	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1179	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
1180	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003919 <+14617>:	movaps 0xb00(%r8),%xmm15
   0x0000000000003921 <+14625>:	mulps  0x380(%rcx),%xmm15
   0x00000000000039df <+14815>:	movaps 0xb40(%r8),%xmm12
   0x00000000000039e7 <+14823>:	mulps  0x380(%rcx),%xmm12
   0x0000000000003a17 <+14871>:	addps  %xmm10,%xmm15
   0x0000000000003a57 <+14935>:	movaps 0xb80(%r8),%xmm14
   0x0000000000003a5f <+14943>:	mulps  0x380(%rcx),%xmm14
   0x0000000000003a7f <+14975>:	addps  %xmm8,%xmm12
   0x0000000000003aab <+15019>:	movaps 0xbc0(%r8),%xmm12
   0x0000000000003ab3 <+15027>:	mulps  0x380(%rcx),%xmm12
   0x0000000000003acf <+15055>:	addps  %xmm7,%xmm14
   0x0000000000003ae6 <+15078>:	addps  %xmm9,%xmm12

1181
1182	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1183	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
1184	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003a1b <+14875>:	movaps 0xb10(%r8),%xmm10
   0x0000000000003a23 <+14883>:	mulps  0x390(%rcx),%xmm10
   0x0000000000003a67 <+14951>:	addps  %xmm15,%xmm10
   0x0000000000003a83 <+14979>:	movaps 0xb50(%r8),%xmm8
   0x0000000000003a8b <+14987>:	mulps  0x390(%rcx),%xmm8
   0x0000000000003aa7 <+15015>:	addps  %xmm12,%xmm8
   0x0000000000003ad3 <+15059>:	movaps 0xb90(%r8),%xmm7
   0x0000000000003adb <+15067>:	mulps  0x390(%rcx),%xmm7
   0x0000000000003aea <+15082>:	movaps 0xbd0(%r8),%xmm9
   0x0000000000003af2 <+15090>:	mulps  0x390(%rcx),%xmm9
   0x0000000000003afa <+15098>:	addps  %xmm14,%xmm7
   0x0000000000003b11 <+15121>:	addps  %xmm12,%xmm9

1185
1186	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1187	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
1188	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003a2f <+14895>:	movaps 0xb60(%r8),%xmm13
   0x0000000000003a37 <+14903>:	mulps  0x3a0(%rcx),%xmm13
   0x0000000000003a43 <+14915>:	movaps 0xb20(%r8),%xmm11
   0x0000000000003a4b <+14923>:	mulps  0x3a0(%rcx),%xmm11
   0x0000000000003a6b <+14955>:	addps  %xmm10,%xmm11
   0x0000000000003a97 <+14999>:	movaps 0xba0(%r8),%xmm11
   0x0000000000003a9f <+15007>:	mulps  0x3a0(%rcx),%xmm11
   0x0000000000003abb <+15035>:	addps  %xmm8,%xmm13
   0x0000000000003afe <+15102>:	addps  %xmm7,%xmm11
   0x0000000000003b29 <+15145>:	movaps 0xbe0(%r8),%xmm11
   0x0000000000003b31 <+15153>:	mulps  0x3a0(%rcx),%xmm11
   0x0000000000003b39 <+15161>:	addps  %xmm9,%xmm11

1189
1190	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1191	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
1192	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003a6f <+14959>:	movaps 0xb30(%r8),%xmm10
   0x0000000000003a77 <+14967>:	mulps  0x3b0(%rcx),%xmm10
   0x0000000000003a93 <+14995>:	addps  %xmm11,%xmm10
   0x0000000000003abf <+15039>:	movaps 0xb70(%r8),%xmm8
   0x0000000000003ac7 <+15047>:	mulps  0x3b0(%rcx),%xmm8
   0x0000000000003ae2 <+15074>:	addps  %xmm13,%xmm8
   0x0000000000003b02 <+15106>:	movaps 0xbb0(%r8),%xmm7
   0x0000000000003b0a <+15114>:	mulps  0x3b0(%rcx),%xmm7
   0x0000000000003b15 <+15125>:	movaps 0xbf0(%r8),%xmm12
   0x0000000000003b1d <+15133>:	mulps  0x3b0(%rcx),%xmm12
   0x0000000000003b25 <+15141>:	addps  %xmm11,%xmm7
   0x0000000000003b3d <+15165>:	addps  %xmm11,%xmm12

1193	      }
1194	    }
1195
1196	    /* Store C(3,3) block. */
1197	    for (i = 0; i < 4; i++)
1198	    {
1199	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000003b41 <+15169>:	mulps  %xmm6,%xmm10
   0x0000000000003b5d <+15197>:	mulps  %xmm6,%xmm8
   0x0000000000003b75 <+15221>:	mulps  %xmm6,%xmm7
   0x0000000000003b8b <+15243>:	mulps  %xmm6,%xmm12

1200	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_33]), C_row[i]);
   0x0000000000003b45 <+15173>:	addps  0x280(%r9),%xmm10
   0x0000000000003b61 <+15201>:	addps  0x290(%r9),%xmm8
   0x0000000000003b78 <+15224>:	addps  0x2a0(%r9),%xmm7
   0x0000000000003b8f <+15247>:	addps  0x2b0(%r9),%xmm12

1201	      _mm_store_ps(&C[i*4+C_OFFSET_33], C_row[i]);
   0x0000000000003b55 <+15189>:	movaps %xmm10,0x280(%r9)
   0x0000000000003b69 <+15209>:	movaps %xmm8,0x290(%r9)
   0x0000000000003b80 <+15232>:	movaps %xmm7,0x2a0(%r9)
   0x0000000000003b97 <+15255>:	movaps %xmm12,0x2b0(%r9)

1202	    }
1203
1204	    /* Reset C(3,4) matrix accumulators */
1205	    C_row[0] = _mm_setzero_ps();
   0x0000000000003b4d <+15181>:	xorps  %xmm9,%xmm9

1206	    C_row[1] = _mm_setzero_ps();
   0x0000000000003b71 <+15217>:	xorps  %xmm8,%xmm8

1207	    C_row[2] = _mm_setzero_ps();
   0x0000000000003b88 <+15240>:	xorps  %xmm7,%xmm7

1208	    C_row[3] = _mm_setzero_ps();
   0x0000000000003b51 <+15185>:	xorps  %xmm11,%xmm11

1209
1210	    if (norm[8]*norm[19] >= tolerance &&
   0x0000000000003b9f <+15263>:	mulss  %xmm2,%xmm0
   0x0000000000003ba3 <+15267>:	comiss %xmm1,%xmm0
   0x0000000000003ba6 <+15270>:	jb     0x40cb <stream_kernel+16587>

1211	        norm[9]*norm[23] >= tolerance &&
   0x0000000000003bac <+15276>:	movss  0x3c(%rdx,%rsi,1),%xmm0
   0x0000000000003bb2 <+15282>:	mulss  0x74(%rdx,%rsi,1),%xmm0
   0x0000000000003bb8 <+15288>:	comiss %xmm1,%xmm0
   0x0000000000003bbb <+15291>:	jb     0x40cb <stream_kernel+16587>

1212	        norm[10]*norm[27] >= tolerance &&
   0x0000000000003bc1 <+15297>:	movss  0x40(%rdx,%rsi,1),%xmm0
   0x0000000000003bc7 <+15303>:	mulss  0x84(%rdx,%rsi,1),%xmm0
   0x0000000000003bd0 <+15312>:	comiss %xmm1,%xmm0
   0x0000000000003bd3 <+15315>:	jb     0x40cb <stream_kernel+16587>

1213	        norm[11]*norm[31] >= tolerance)
   0x0000000000003bd9 <+15321>:	movss  0x44(%rdx,%rsi,1),%xmm0
   0x0000000000003bdf <+15327>:	mulss  0x94(%rdx,%rsi,1),%xmm0
   0x0000000000003be8 <+15336>:	comiss %xmm1,%xmm0
   0x0000000000003beb <+15339>:	jb     0x40cb <stream_kernel+16587>

1214	    {
1215	      /* A(3,1)*B(1,4) = C(3,4). */
1216	      for (i = 0; i < 4; i++)
1217	      {
1218	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_31]);
1219	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
1220	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003bf1 <+15345>:	movaps 0x800(%r8),%xmm7
   0x0000000000003c11 <+15377>:	movaps 0x840(%r8),%xmm11
   0x0000000000003c31 <+15409>:	movaps 0x880(%r8),%xmm14
   0x0000000000003c41 <+15425>:	mulps  0xc0(%rcx),%xmm7
   0x0000000000003c5f <+15455>:	mulps  0xc0(%rcx),%xmm11
   0x0000000000003c7f <+15487>:	mulps  0xc0(%rcx),%xmm14
   0x0000000000003cb8 <+15544>:	movaps 0x8c0(%r8),%xmm10
   0x0000000000003cc0 <+15552>:	mulps  0xc0(%rcx),%xmm10

1221
1222	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_31]);
1223	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
1224	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003bf9 <+15353>:	movaps 0x810(%r8),%xmm0
   0x0000000000003c19 <+15385>:	movaps 0x850(%r8),%xmm12
   0x0000000000003c48 <+15432>:	mulps  0xd0(%rcx),%xmm0
   0x0000000000003c67 <+15463>:	mulps  0xd0(%rcx),%xmm12
   0x0000000000003c8f <+15503>:	addps  %xmm7,%xmm0
   0x0000000000003c92 <+15506>:	movaps 0x890(%r8),%xmm7
   0x0000000000003c9a <+15514>:	mulps  0xd0(%rcx),%xmm7
   0x0000000000003cc8 <+15560>:	addps  %xmm11,%xmm12
   0x0000000000003ccc <+15564>:	movaps 0x8d0(%r8),%xmm11
   0x0000000000003cd4 <+15572>:	mulps  0xd0(%rcx),%xmm11
   0x0000000000003d04 <+15620>:	addps  %xmm14,%xmm7
   0x0000000000003d3d <+15677>:	addps  %xmm10,%xmm11

1225
1226	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_31]);
1227	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
1228	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003c01 <+15361>:	movaps 0x820(%r8),%xmm10
   0x0000000000003c21 <+15393>:	movaps 0x860(%r8),%xmm13
   0x0000000000003c4f <+15439>:	mulps  0xe0(%rcx),%xmm10
   0x0000000000003c6f <+15471>:	mulps  0xe0(%rcx),%xmm13
   0x0000000000003ca1 <+15521>:	addps  %xmm0,%xmm10
   0x0000000000003ca5 <+15525>:	movaps 0x8a0(%r8),%xmm0
   0x0000000000003cad <+15533>:	mulps  0xe0(%rcx),%xmm0
   0x0000000000003cdc <+15580>:	addps  %xmm12,%xmm13
   0x0000000000003ce0 <+15584>:	movaps 0x8e0(%r8),%xmm12
   0x0000000000003ce8 <+15592>:	mulps  0xe0(%rcx),%xmm12
   0x0000000000003d18 <+15640>:	addps  %xmm7,%xmm0
   0x0000000000003d51 <+15697>:	addps  %xmm11,%xmm12

1229
1230	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_31]);
1231	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
1232	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003c09 <+15369>:	movaps 0x830(%r8),%xmm8
   0x0000000000003c29 <+15401>:	movaps 0x870(%r8),%xmm9
   0x0000000000003c39 <+15417>:	movaps 0x8b0(%r8),%xmm15
   0x0000000000003c57 <+15447>:	mulps  0xf0(%rcx),%xmm8
   0x0000000000003c77 <+15479>:	mulps  0xf0(%rcx),%xmm9
   0x0000000000003c87 <+15495>:	mulps  0xf0(%rcx),%xmm15
   0x0000000000003cb4 <+15540>:	addps  %xmm10,%xmm8
   0x0000000000003cf0 <+15600>:	addps  %xmm13,%xmm9
   0x0000000000003d08 <+15624>:	movaps 0x8f0(%r8),%xmm14
   0x0000000000003d10 <+15632>:	mulps  0xf0(%rcx),%xmm14
   0x0000000000003d2a <+15658>:	addps  %xmm0,%xmm15
   0x0000000000003d65 <+15717>:	addps  %xmm12,%xmm14

1233	      }
1234
1235	      /* A(3,2)*B(2,4) = C(3,4). */
1236	      for (i = 0; i < 4; i++)
1237	      {
1238	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_32]);
1239	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
1240	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003d1b <+15643>:	movaps 0x900(%r8),%xmm7
   0x0000000000003d23 <+15651>:	mulps  0x1c0(%rcx),%xmm7
   0x0000000000003d2e <+15662>:	movaps 0x940(%r8),%xmm0
   0x0000000000003d36 <+15670>:	mulps  0x1c0(%rcx),%xmm0
   0x0000000000003d79 <+15737>:	addps  %xmm8,%xmm7
   0x0000000000003d8d <+15757>:	addps  %xmm9,%xmm0
   0x0000000000003ddf <+15839>:	movaps 0x980(%r8),%xmm7
   0x0000000000003de7 <+15847>:	mulps  0x1c0(%rcx),%xmm7
   0x0000000000003e06 <+15878>:	movaps 0x9c0(%r8),%xmm8
   0x0000000000003e0e <+15886>:	mulps  0x1c0(%rcx),%xmm8
   0x0000000000003e16 <+15894>:	addps  %xmm15,%xmm7
   0x0000000000003e2a <+15914>:	addps  %xmm14,%xmm8

1241
1242	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_32]);
1243	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
1244	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003d7d <+15741>:	movaps 0x910(%r8),%xmm8
   0x0000000000003d85 <+15749>:	mulps  0x1d0(%rcx),%xmm8
   0x0000000000003d91 <+15761>:	movaps 0x950(%r8),%xmm9
   0x0000000000003d99 <+15769>:	mulps  0x1d0(%rcx),%xmm9
   0x0000000000003da1 <+15777>:	addps  %xmm7,%xmm8
   0x0000000000003db4 <+15796>:	addps  %xmm0,%xmm9
   0x0000000000003e1a <+15898>:	movaps 0x990(%r8),%xmm15
   0x0000000000003e22 <+15906>:	mulps  0x1d0(%rcx),%xmm15
   0x0000000000003e3e <+15934>:	addps  %xmm7,%xmm15
   0x0000000000003e42 <+15938>:	movaps 0x9d0(%r8),%xmm7
   0x0000000000003e4a <+15946>:	mulps  0x1d0(%rcx),%xmm7
   0x0000000000003e79 <+15993>:	addps  %xmm8,%xmm7

1245
1246	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_32]);
1247	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
1248	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003da5 <+15781>:	movaps 0x920(%r8),%xmm7
   0x0000000000003dad <+15789>:	mulps  0x1e0(%rcx),%xmm7
   0x0000000000003db8 <+15800>:	movaps 0x9e0(%r8),%xmm0
   0x0000000000003dc0 <+15808>:	mulps  0x1e0(%rcx),%xmm0
   0x0000000000003dc7 <+15815>:	addps  %xmm8,%xmm7
   0x0000000000003dcb <+15819>:	movaps 0x960(%r8),%xmm8
   0x0000000000003dd3 <+15827>:	mulps  0x1e0(%rcx),%xmm8
   0x0000000000003dee <+15854>:	addps  %xmm9,%xmm8
   0x0000000000003df2 <+15858>:	movaps 0x9a0(%r8),%xmm9
   0x0000000000003dfa <+15866>:	mulps  0x1e0(%rcx),%xmm9
   0x0000000000003e51 <+15953>:	addps  %xmm15,%xmm9
   0x0000000000003e8d <+16013>:	addps  %xmm7,%xmm0

1249
1250	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_32]);
1251	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
1252	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003cf4 <+15604>:	movaps 0x930(%r8),%xmm13
   0x0000000000003cfc <+15612>:	mulps  0x1f0(%rcx),%xmm13
   0x0000000000003d41 <+15681>:	movaps 0x9f0(%r8),%xmm10
   0x0000000000003d49 <+15689>:	mulps  0x1f0(%rcx),%xmm10
   0x0000000000003d55 <+15701>:	movaps 0x9b0(%r8),%xmm11
   0x0000000000003d5d <+15709>:	mulps  0x1f0(%rcx),%xmm11
   0x0000000000003d69 <+15721>:	movaps 0x970(%r8),%xmm12
   0x0000000000003d71 <+15729>:	mulps  0x1f0(%rcx),%xmm12
   0x0000000000003ddb <+15835>:	addps  %xmm7,%xmm13
   0x0000000000003e02 <+15874>:	addps  %xmm8,%xmm12
   0x0000000000003e65 <+15973>:	addps  %xmm9,%xmm11
   0x0000000000003e9f <+16031>:	addps  %xmm0,%xmm10

1253	      }
1254
1255	      /* A(3,3)*B(3,4) = C(3,4). */
1256	      for (i = 0; i < 4; i++)
1257	      {
1258	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_33]);
1259	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
1260	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003e2e <+15918>:	movaps 0xa00(%r8),%xmm14
   0x0000000000003e36 <+15926>:	mulps  0x2c0(%rcx),%xmm14
   0x0000000000003e90 <+16016>:	movaps 0xa40(%r8),%xmm7
   0x0000000000003e98 <+16024>:	mulps  0x2c0(%rcx),%xmm7
   0x0000000000003eb2 <+16050>:	addps  %xmm13,%xmm14
   0x0000000000003ec6 <+16070>:	addps  %xmm12,%xmm7
   0x0000000000003ede <+16094>:	movaps 0xa80(%r8),%xmm14
   0x0000000000003ee6 <+16102>:	mulps  0x2c0(%rcx),%xmm14
   0x0000000000003f2c <+16172>:	movaps 0xac0(%r8),%xmm0
   0x0000000000003f34 <+16180>:	mulps  0x2c0(%rcx),%xmm0
   0x0000000000003f4f <+16207>:	addps  %xmm11,%xmm14
   0x0000000000003f63 <+16227>:	addps  %xmm10,%xmm0

1261
1262	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_33]);
1263	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
1264	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003ea3 <+16035>:	movaps 0xa50(%r8),%xmm0
   0x0000000000003eab <+16043>:	mulps  0x2d0(%rcx),%xmm0
   0x0000000000003eb6 <+16054>:	movaps 0xa10(%r8),%xmm13
   0x0000000000003ebe <+16062>:	mulps  0x2d0(%rcx),%xmm13
   0x0000000000003eda <+16090>:	addps  %xmm14,%xmm13
   0x0000000000003f16 <+16150>:	addps  %xmm7,%xmm0
   0x0000000000003f53 <+16211>:	movaps 0xa90(%r8),%xmm11
   0x0000000000003f5b <+16219>:	mulps  0x2d0(%rcx),%xmm11
   0x0000000000003f67 <+16231>:	movaps 0xad0(%r8),%xmm10
   0x0000000000003f6f <+16239>:	mulps  0x2d0(%rcx),%xmm10
   0x0000000000003f77 <+16247>:	addps  %xmm14,%xmm11
   0x0000000000003fb3 <+16307>:	addps  %xmm0,%xmm10

1265
1266	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_33]);
1267	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
1268	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003e55 <+15957>:	movaps 0xaa0(%r8),%xmm15
   0x0000000000003e5d <+15965>:	mulps  0x2e0(%rcx),%xmm15
   0x0000000000003e7d <+15997>:	movaps 0xa20(%r8),%xmm8
   0x0000000000003e85 <+16005>:	mulps  0x2e0(%rcx),%xmm8
   0x0000000000003eca <+16074>:	movaps 0xa60(%r8),%xmm12
   0x0000000000003ed2 <+16082>:	mulps  0x2e0(%rcx),%xmm12
   0x0000000000003eee <+16110>:	addps  %xmm13,%xmm8
   0x0000000000003f28 <+16168>:	addps  %xmm0,%xmm12
   0x0000000000003f8b <+16267>:	addps  %xmm11,%xmm15
   0x0000000000003f8f <+16271>:	movaps 0xae0(%r8),%xmm11
   0x0000000000003f97 <+16279>:	mulps  0x2e0(%rcx),%xmm11
   0x0000000000003fc6 <+16326>:	addps  %xmm10,%xmm11

1269
1270	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_33]);
1271	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
1272	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003e69 <+15977>:	movaps 0xa30(%r8),%xmm9
   0x0000000000003e71 <+15985>:	mulps  0x2f0(%rcx),%xmm9
   0x0000000000003f02 <+16130>:	addps  %xmm8,%xmm9
   0x0000000000003f06 <+16134>:	movaps 0xa70(%r8),%xmm8
   0x0000000000003f0e <+16142>:	mulps  0x2f0(%rcx),%xmm8
   0x0000000000003f19 <+16153>:	movaps 0xab0(%r8),%xmm7
   0x0000000000003f21 <+16161>:	mulps  0x2f0(%rcx),%xmm7
   0x0000000000003f3b <+16187>:	addps  %xmm12,%xmm8
   0x0000000000003f9f <+16287>:	addps  %xmm15,%xmm7
   0x0000000000003fb7 <+16311>:	movaps 0xaf0(%r8),%xmm0
   0x0000000000003fbf <+16319>:	mulps  0x2f0(%rcx),%xmm0
   0x0000000000003fda <+16346>:	addps  %xmm11,%xmm0

1273	      }
1274
1275	      /* A(3,4)*B(4,4) = C(3,4). */
1276	      for (i = 0; i < 4; i++)
1277	      {
1278	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_34]);
1279	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
1280	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003f3f <+16191>:	movaps 0xb00(%r8),%xmm12
   0x0000000000003f47 <+16199>:	mulps  0x3c0(%rcx),%xmm12
   0x0000000000003f7b <+16251>:	movaps 0xb40(%r8),%xmm14
   0x0000000000003f83 <+16259>:	mulps  0x3c0(%rcx),%xmm14
   0x0000000000003fca <+16330>:	movaps 0xb80(%r8),%xmm10
   0x0000000000003fd2 <+16338>:	mulps  0x3c0(%rcx),%xmm10
   0x0000000000003fee <+16366>:	addps  %xmm9,%xmm12
   0x0000000000004002 <+16386>:	addps  %xmm8,%xmm14
   0x000000000000401a <+16410>:	movaps 0xbc0(%r8),%xmm12
   0x0000000000004022 <+16418>:	mulps  0x3c0(%rcx),%xmm12
   0x000000000000405a <+16474>:	addps  %xmm7,%xmm10
   0x0000000000004071 <+16497>:	addps  %xmm0,%xmm12

1281
1282	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_34]);
1283	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
1284	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003ff2 <+16370>:	movaps 0xb10(%r8),%xmm9
   0x0000000000003ffa <+16378>:	mulps  0x3d0(%rcx),%xmm9
   0x0000000000004006 <+16390>:	movaps 0xb50(%r8),%xmm8
   0x000000000000400e <+16398>:	mulps  0x3d0(%rcx),%xmm8
   0x0000000000004016 <+16406>:	addps  %xmm12,%xmm9
   0x000000000000403e <+16446>:	addps  %xmm14,%xmm8
   0x000000000000405e <+16478>:	movaps 0xb90(%r8),%xmm7
   0x0000000000004066 <+16486>:	mulps  0x3d0(%rcx),%xmm7
   0x0000000000004075 <+16501>:	movaps 0xbd0(%r8),%xmm0
   0x000000000000407d <+16509>:	mulps  0x3d0(%rcx),%xmm0
   0x0000000000004084 <+16516>:	addps  %xmm10,%xmm7
   0x00000000000040ab <+16555>:	addps  %xmm12,%xmm0

1285
1286	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_34]);
1287	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
1288	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000003ef2 <+16114>:	movaps 0xb20(%r8),%xmm13
   0x0000000000003efa <+16122>:	mulps  0x3e0(%rcx),%xmm13
   0x0000000000003fa3 <+16291>:	movaps 0xb60(%r8),%xmm15
   0x0000000000003fab <+16299>:	mulps  0x3e0(%rcx),%xmm15
   0x0000000000003fde <+16350>:	movaps 0xba0(%r8),%xmm11
   0x0000000000003fe6 <+16358>:	mulps  0x3e0(%rcx),%xmm11
   0x000000000000402a <+16426>:	addps  %xmm9,%xmm13
   0x0000000000004042 <+16450>:	addps  %xmm8,%xmm15
   0x0000000000004088 <+16520>:	movaps 0xbe0(%r8),%xmm10
   0x0000000000004090 <+16528>:	mulps  0x3e0(%rcx),%xmm10
   0x0000000000004098 <+16536>:	addps  %xmm7,%xmm11
   0x00000000000040af <+16559>:	addps  %xmm0,%xmm10

1289
1290	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_34]);
1291	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
1292	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000402e <+16430>:	movaps 0xb30(%r8),%xmm9
   0x0000000000004036 <+16438>:	mulps  0x3f0(%rcx),%xmm9
   0x0000000000004046 <+16454>:	movaps 0xb70(%r8),%xmm8
   0x000000000000404e <+16462>:	mulps  0x3f0(%rcx),%xmm8
   0x0000000000004056 <+16470>:	addps  %xmm13,%xmm9
   0x000000000000406d <+16493>:	addps  %xmm15,%xmm8
   0x000000000000409c <+16540>:	movaps 0xbb0(%r8),%xmm7
   0x00000000000040a4 <+16548>:	mulps  0x3f0(%rcx),%xmm7
   0x00000000000040b3 <+16563>:	addps  %xmm11,%xmm7
   0x00000000000040b7 <+16567>:	movaps 0xbf0(%r8),%xmm11
   0x00000000000040bf <+16575>:	mulps  0x3f0(%rcx),%xmm11
   0x00000000000040c7 <+16583>:	addps  %xmm10,%xmm11

1293	      }
1294	    }
1295
1296	    /* Store C(3,4) block. */
1297	    for (i = 0; i < 4; i++)
1298	    {
1299	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x00000000000040cb <+16587>:	mulps  %xmm6,%xmm9
   0x00000000000040e9 <+16617>:	mulps  %xmm6,%xmm8
   0x0000000000004101 <+16641>:	mulps  %xmm6,%xmm7
   0x0000000000004117 <+16663>:	mulps  %xmm6,%xmm11

1300	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_34]), C_row[i]);
   0x00000000000040cf <+16591>:	addps  0x2c0(%r9),%xmm9
   0x00000000000040ed <+16621>:	addps  0x2d0(%r9),%xmm8
   0x0000000000004104 <+16644>:	addps  0x2e0(%r9),%xmm7
   0x000000000000411b <+16667>:	addps  0x2f0(%r9),%xmm11

1301	      _mm_store_ps(&C[i*4+C_OFFSET_34], C_row[i]);
   0x00000000000040dd <+16605>:	movaps %xmm9,0x2c0(%r9)
   0x00000000000040f5 <+16629>:	movaps %xmm8,0x2d0(%r9)
   0x000000000000410c <+16652>:	movaps %xmm7,0x2e0(%r9)
   0x0000000000004123 <+16675>:	movaps %xmm11,0x2f0(%r9)

1302	    }
1303
1304	    /* Reset C(4,1) matrix accumulators */
1305	    C_row[0] = _mm_setzero_ps();
   0x00000000000040e5 <+16613>:	xorps  %xmm9,%xmm9

1306	    C_row[1] = _mm_setzero_ps();
   0x00000000000040fd <+16637>:	xorps  %xmm8,%xmm8

1307	    C_row[2] = _mm_setzero_ps();
   0x0000000000004114 <+16660>:	xorps  %xmm7,%xmm7

1308	    C_row[3] = _mm_setzero_ps();
   0x000000000000412b <+16683>:	xorps  %xmm11,%xmm11

1309
1310	    if (norm[12]*norm[16] >= tolerance &&
   0x00000000000040d7 <+16599>:	movss  0x48(%rdx,%rsi,1),%xmm0
   0x000000000000412f <+16687>:	mulss  %xmm0,%xmm5
   0x0000000000004133 <+16691>:	comiss %xmm1,%xmm5
   0x0000000000004136 <+16694>:	jb     0x4624 <stream_kernel+17956>

1311	        norm[13]*norm[20] >= tolerance &&
   0x000000000000413c <+16700>:	movss  0x4c(%rdx,%rsi,1),%xmm5
   0x0000000000004142 <+16706>:	mulss  0x68(%rdx,%rsi,1),%xmm5
   0x0000000000004148 <+16712>:	comiss %xmm1,%xmm5
   0x000000000000414b <+16715>:	jb     0x4624 <stream_kernel+17956>

1312	        norm[14]*norm[24] >= tolerance &&
   0x0000000000004151 <+16721>:	movss  0x50(%rdx,%rsi,1),%xmm5
   0x0000000000004157 <+16727>:	mulss  0x78(%rdx,%rsi,1),%xmm5
   0x000000000000415d <+16733>:	comiss %xmm1,%xmm5
   0x0000000000004160 <+16736>:	jb     0x4624 <stream_kernel+17956>

1313	        norm[15]*norm[28] >= tolerance)
   0x0000000000004166 <+16742>:	movss  0x54(%rdx,%rsi,1),%xmm5
   0x000000000000416c <+16748>:	mulss  0x88(%rdx,%rsi,1),%xmm5
   0x0000000000004175 <+16757>:	comiss %xmm1,%xmm5
   0x0000000000004178 <+16760>:	jb     0x4624 <stream_kernel+17956>

1314	    {
1315	      /* A(4,1)*B(1,1) = C(4,1). */
1316	      for (i = 0; i < 4; i++)
1317	      {
1318	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1319	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
1320	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000417e <+16766>:	movaps 0xc00(%r8),%xmm7
   0x000000000000419e <+16798>:	movaps 0xc40(%r8),%xmm11
   0x00000000000041be <+16830>:	movaps 0xc80(%r8),%xmm14
   0x00000000000041ce <+16846>:	mulps  (%rcx),%xmm7
   0x00000000000041df <+16863>:	mulps  (%rcx),%xmm11
   0x00000000000041f2 <+16882>:	mulps  (%rcx),%xmm14
   0x000000000000421e <+16926>:	movaps 0xcc0(%r8),%xmm10
   0x0000000000004226 <+16934>:	mulps  (%rcx),%xmm10

1321
1322	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1323	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
1324	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004186 <+16774>:	movaps 0xc10(%r8),%xmm5
   0x00000000000041a6 <+16806>:	movaps 0xc50(%r8),%xmm12
   0x00000000000041d1 <+16849>:	mulps  0x10(%rcx),%xmm5
   0x00000000000041e3 <+16867>:	mulps  0x10(%rcx),%xmm12
   0x00000000000041fb <+16891>:	addps  %xmm7,%xmm5
   0x00000000000041fe <+16894>:	movaps 0xc90(%r8),%xmm7
   0x0000000000004206 <+16902>:	mulps  0x10(%rcx),%xmm7
   0x000000000000422a <+16938>:	addps  %xmm11,%xmm12
   0x000000000000422e <+16942>:	movaps 0xcd0(%r8),%xmm11
   0x0000000000004236 <+16950>:	mulps  0x10(%rcx),%xmm11
   0x0000000000004260 <+16992>:	addps  %xmm14,%xmm7
   0x0000000000004296 <+17046>:	addps  %xmm10,%xmm11

1325
1326	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1327	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
1328	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000418e <+16782>:	movaps 0xc20(%r8),%xmm10
   0x00000000000041ae <+16814>:	movaps 0xc60(%r8),%xmm13
   0x00000000000041d5 <+16853>:	mulps  0x20(%rcx),%xmm10
   0x00000000000041e8 <+16872>:	mulps  0x20(%rcx),%xmm13
   0x000000000000420a <+16906>:	addps  %xmm5,%xmm10
   0x000000000000420e <+16910>:	movaps 0xca0(%r8),%xmm5
   0x0000000000004216 <+16918>:	mulps  0x20(%rcx),%xmm5
   0x000000000000423b <+16955>:	addps  %xmm12,%xmm13
   0x000000000000423f <+16959>:	movaps 0xce0(%r8),%xmm12
   0x0000000000004247 <+16967>:	mulps  0x20(%rcx),%xmm12
   0x0000000000004271 <+17009>:	addps  %xmm7,%xmm5
   0x00000000000042aa <+17066>:	addps  %xmm11,%xmm12

1329
1330	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1331	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
1332	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004196 <+16790>:	movaps 0xc30(%r8),%xmm8
   0x00000000000041b6 <+16822>:	movaps 0xc70(%r8),%xmm9
   0x00000000000041c6 <+16838>:	movaps 0xcb0(%r8),%xmm15
   0x00000000000041da <+16858>:	mulps  0x30(%rcx),%xmm8
   0x00000000000041ed <+16877>:	mulps  0x30(%rcx),%xmm9
   0x00000000000041f6 <+16886>:	mulps  0x30(%rcx),%xmm15
   0x000000000000421a <+16922>:	addps  %xmm10,%xmm8
   0x000000000000424c <+16972>:	addps  %xmm13,%xmm9
   0x0000000000004264 <+16996>:	movaps 0xcf0(%r8),%xmm14
   0x000000000000426c <+17004>:	mulps  0x30(%rcx),%xmm14
   0x0000000000004283 <+17027>:	addps  %xmm5,%xmm15
   0x00000000000042be <+17086>:	addps  %xmm12,%xmm14

1333	      }
1334
1335	      /* A(4,2)*B(2,1) = C(4,1). */
1336	      for (i = 0; i < 4; i++)
1337	      {
1338	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1339	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
1340	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004274 <+17012>:	movaps 0xd00(%r8),%xmm7
   0x000000000000427c <+17020>:	mulps  0x100(%rcx),%xmm7
   0x0000000000004287 <+17031>:	movaps 0xd40(%r8),%xmm5
   0x000000000000428f <+17039>:	mulps  0x100(%rcx),%xmm5
   0x00000000000042d2 <+17106>:	addps  %xmm8,%xmm7
   0x00000000000042e6 <+17126>:	addps  %xmm9,%xmm5
   0x0000000000004338 <+17208>:	movaps 0xd80(%r8),%xmm7
   0x0000000000004340 <+17216>:	mulps  0x100(%rcx),%xmm7
   0x000000000000435f <+17247>:	movaps 0xdc0(%r8),%xmm8
   0x0000000000004367 <+17255>:	mulps  0x100(%rcx),%xmm8
   0x000000000000436f <+17263>:	addps  %xmm15,%xmm7
   0x0000000000004383 <+17283>:	addps  %xmm14,%xmm8

1341
1342	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1343	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
1344	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000042d6 <+17110>:	movaps 0xd10(%r8),%xmm8
   0x00000000000042de <+17118>:	mulps  0x110(%rcx),%xmm8
   0x00000000000042ea <+17130>:	movaps 0xd50(%r8),%xmm9
   0x00000000000042f2 <+17138>:	mulps  0x110(%rcx),%xmm9
   0x00000000000042fa <+17146>:	addps  %xmm7,%xmm8
   0x000000000000430d <+17165>:	addps  %xmm5,%xmm9
   0x0000000000004373 <+17267>:	movaps 0xd90(%r8),%xmm15
   0x000000000000437b <+17275>:	mulps  0x110(%rcx),%xmm15
   0x0000000000004397 <+17303>:	addps  %xmm7,%xmm15
   0x000000000000439b <+17307>:	movaps 0xdd0(%r8),%xmm7
   0x00000000000043a3 <+17315>:	mulps  0x110(%rcx),%xmm7
   0x00000000000043d2 <+17362>:	addps  %xmm8,%xmm7

1345
1346	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1347	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
1348	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000042fe <+17150>:	movaps 0xd20(%r8),%xmm7
   0x0000000000004306 <+17158>:	mulps  0x120(%rcx),%xmm7
   0x0000000000004311 <+17169>:	movaps 0xde0(%r8),%xmm5
   0x0000000000004319 <+17177>:	mulps  0x120(%rcx),%xmm5
   0x0000000000004320 <+17184>:	addps  %xmm8,%xmm7
   0x0000000000004324 <+17188>:	movaps 0xd60(%r8),%xmm8
   0x000000000000432c <+17196>:	mulps  0x120(%rcx),%xmm8
   0x0000000000004347 <+17223>:	addps  %xmm9,%xmm8
   0x000000000000434b <+17227>:	movaps 0xda0(%r8),%xmm9
   0x0000000000004353 <+17235>:	mulps  0x120(%rcx),%xmm9
   0x00000000000043aa <+17322>:	addps  %xmm15,%xmm9
   0x00000000000043e6 <+17382>:	addps  %xmm7,%xmm5

1349
1350	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1351	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
1352	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004250 <+16976>:	movaps 0xd30(%r8),%xmm13
   0x0000000000004258 <+16984>:	mulps  0x130(%rcx),%xmm13
   0x000000000000429a <+17050>:	movaps 0xdf0(%r8),%xmm10
   0x00000000000042a2 <+17058>:	mulps  0x130(%rcx),%xmm10
   0x00000000000042ae <+17070>:	movaps 0xdb0(%r8),%xmm11
   0x00000000000042b6 <+17078>:	mulps  0x130(%rcx),%xmm11
   0x00000000000042c2 <+17090>:	movaps 0xd70(%r8),%xmm12
   0x00000000000042ca <+17098>:	mulps  0x130(%rcx),%xmm12
   0x0000000000004334 <+17204>:	addps  %xmm7,%xmm13
   0x000000000000435b <+17243>:	addps  %xmm8,%xmm12
   0x00000000000043be <+17342>:	addps  %xmm9,%xmm11
   0x00000000000043f8 <+17400>:	addps  %xmm5,%xmm10

1353	      }
1354
1355	      /* A(4,3)*B(3,1) = C(4,1). */
1356	      for (i = 0; i < 4; i++)
1357	      {
1358	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1359	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_31]);
1360	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004387 <+17287>:	movaps 0xe00(%r8),%xmm14
   0x000000000000438f <+17295>:	mulps  0x200(%rcx),%xmm14
   0x00000000000043e9 <+17385>:	movaps 0xe40(%r8),%xmm7
   0x00000000000043f1 <+17393>:	mulps  0x200(%rcx),%xmm7
   0x000000000000440b <+17419>:	addps  %xmm13,%xmm14
   0x000000000000441f <+17439>:	addps  %xmm12,%xmm7
   0x0000000000004437 <+17463>:	movaps 0xe80(%r8),%xmm14
   0x000000000000443f <+17471>:	mulps  0x200(%rcx),%xmm14
   0x0000000000004485 <+17541>:	movaps 0xec0(%r8),%xmm5
   0x000000000000448d <+17549>:	mulps  0x200(%rcx),%xmm5
   0x00000000000044a8 <+17576>:	addps  %xmm11,%xmm14
   0x00000000000044bc <+17596>:	addps  %xmm10,%xmm5

1361
1362	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1363	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_31]);
1364	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000043fc <+17404>:	movaps 0xe50(%r8),%xmm5
   0x0000000000004404 <+17412>:	mulps  0x210(%rcx),%xmm5
   0x000000000000440f <+17423>:	movaps 0xe10(%r8),%xmm13
   0x0000000000004417 <+17431>:	mulps  0x210(%rcx),%xmm13
   0x0000000000004433 <+17459>:	addps  %xmm14,%xmm13
   0x000000000000446f <+17519>:	addps  %xmm7,%xmm5
   0x00000000000044ac <+17580>:	movaps 0xe90(%r8),%xmm11
   0x00000000000044b4 <+17588>:	mulps  0x210(%rcx),%xmm11
   0x00000000000044c0 <+17600>:	movaps 0xed0(%r8),%xmm10
   0x00000000000044c8 <+17608>:	mulps  0x210(%rcx),%xmm10
   0x00000000000044d0 <+17616>:	addps  %xmm14,%xmm11
   0x000000000000450c <+17676>:	addps  %xmm5,%xmm10

1365
1366	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1367	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_31]);
1368	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000043ae <+17326>:	movaps 0xea0(%r8),%xmm15
   0x00000000000043b6 <+17334>:	mulps  0x220(%rcx),%xmm15
   0x00000000000043d6 <+17366>:	movaps 0xe20(%r8),%xmm8
   0x00000000000043de <+17374>:	mulps  0x220(%rcx),%xmm8
   0x0000000000004423 <+17443>:	movaps 0xe60(%r8),%xmm12
   0x000000000000442b <+17451>:	mulps  0x220(%rcx),%xmm12
   0x0000000000004447 <+17479>:	addps  %xmm13,%xmm8
   0x0000000000004481 <+17537>:	addps  %xmm5,%xmm12
   0x00000000000044e4 <+17636>:	addps  %xmm11,%xmm15
   0x00000000000044e8 <+17640>:	movaps 0xee0(%r8),%xmm11
   0x00000000000044f0 <+17648>:	mulps  0x220(%rcx),%xmm11
   0x000000000000451f <+17695>:	addps  %xmm10,%xmm11

1369
1370	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1371	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_31]);
1372	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000043c2 <+17346>:	movaps 0xe30(%r8),%xmm9
   0x00000000000043ca <+17354>:	mulps  0x230(%rcx),%xmm9
   0x000000000000445b <+17499>:	addps  %xmm8,%xmm9
   0x000000000000445f <+17503>:	movaps 0xe70(%r8),%xmm8
   0x0000000000004467 <+17511>:	mulps  0x230(%rcx),%xmm8
   0x0000000000004472 <+17522>:	movaps 0xeb0(%r8),%xmm7
   0x000000000000447a <+17530>:	mulps  0x230(%rcx),%xmm7
   0x0000000000004494 <+17556>:	addps  %xmm12,%xmm8
   0x00000000000044f8 <+17656>:	addps  %xmm15,%xmm7
   0x0000000000004510 <+17680>:	movaps 0xef0(%r8),%xmm5
   0x0000000000004518 <+17688>:	mulps  0x230(%rcx),%xmm5
   0x0000000000004533 <+17715>:	addps  %xmm11,%xmm5

1373	      }
1374
1375	      /* A(4,4)*B(4,1) = C(4,1). */
1376	      for (i = 0; i < 4; i++)
1377	      {
1378	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1379	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_41]);
1380	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004498 <+17560>:	movaps 0xf00(%r8),%xmm12
   0x00000000000044a0 <+17568>:	mulps  0x300(%rcx),%xmm12
   0x00000000000044d4 <+17620>:	movaps 0xf40(%r8),%xmm14
   0x00000000000044dc <+17628>:	mulps  0x300(%rcx),%xmm14
   0x0000000000004523 <+17699>:	movaps 0xf80(%r8),%xmm10
   0x000000000000452b <+17707>:	mulps  0x300(%rcx),%xmm10
   0x0000000000004547 <+17735>:	addps  %xmm9,%xmm12
   0x000000000000455b <+17755>:	addps  %xmm8,%xmm14
   0x0000000000004573 <+17779>:	movaps 0xfc0(%r8),%xmm12
   0x000000000000457b <+17787>:	mulps  0x300(%rcx),%xmm12
   0x00000000000045b3 <+17843>:	addps  %xmm7,%xmm10
   0x00000000000045ca <+17866>:	addps  %xmm5,%xmm12

1381
1382	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1383	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_41]);
1384	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000454b <+17739>:	movaps 0xf10(%r8),%xmm9
   0x0000000000004553 <+17747>:	mulps  0x310(%rcx),%xmm9
   0x000000000000455f <+17759>:	movaps 0xf50(%r8),%xmm8
   0x0000000000004567 <+17767>:	mulps  0x310(%rcx),%xmm8
   0x000000000000456f <+17775>:	addps  %xmm12,%xmm9
   0x0000000000004597 <+17815>:	addps  %xmm14,%xmm8
   0x00000000000045b7 <+17847>:	movaps 0xf90(%r8),%xmm7
   0x00000000000045bf <+17855>:	mulps  0x310(%rcx),%xmm7
   0x00000000000045ce <+17870>:	movaps 0xfd0(%r8),%xmm5
   0x00000000000045d6 <+17878>:	mulps  0x310(%rcx),%xmm5
   0x00000000000045dd <+17885>:	addps  %xmm10,%xmm7
   0x0000000000004604 <+17924>:	addps  %xmm12,%xmm5

1385
1386	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1387	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_41]);
1388	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000444b <+17483>:	movaps 0xf20(%r8),%xmm13
   0x0000000000004453 <+17491>:	mulps  0x320(%rcx),%xmm13
   0x00000000000044fc <+17660>:	movaps 0xf60(%r8),%xmm15
   0x0000000000004504 <+17668>:	mulps  0x320(%rcx),%xmm15
   0x0000000000004537 <+17719>:	movaps 0xfa0(%r8),%xmm11
   0x000000000000453f <+17727>:	mulps  0x320(%rcx),%xmm11
   0x0000000000004583 <+17795>:	addps  %xmm9,%xmm13
   0x000000000000459b <+17819>:	addps  %xmm8,%xmm15
   0x00000000000045e1 <+17889>:	movaps 0xfe0(%r8),%xmm10
   0x00000000000045e9 <+17897>:	mulps  0x320(%rcx),%xmm10
   0x00000000000045f1 <+17905>:	addps  %xmm7,%xmm11
   0x0000000000004608 <+17928>:	addps  %xmm5,%xmm10

1389
1390	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1391	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_41]);
1392	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004587 <+17799>:	movaps 0xf30(%r8),%xmm9
   0x000000000000458f <+17807>:	mulps  0x330(%rcx),%xmm9
   0x000000000000459f <+17823>:	movaps 0xf70(%r8),%xmm8
   0x00000000000045a7 <+17831>:	mulps  0x330(%rcx),%xmm8
   0x00000000000045af <+17839>:	addps  %xmm13,%xmm9
   0x00000000000045c6 <+17862>:	addps  %xmm15,%xmm8
   0x00000000000045f5 <+17909>:	movaps 0xfb0(%r8),%xmm7
   0x00000000000045fd <+17917>:	mulps  0x330(%rcx),%xmm7
   0x000000000000460c <+17932>:	addps  %xmm11,%xmm7
   0x0000000000004610 <+17936>:	movaps 0xff0(%r8),%xmm11
   0x0000000000004618 <+17944>:	mulps  0x330(%rcx),%xmm11
   0x0000000000004620 <+17952>:	addps  %xmm10,%xmm11

1393	      }
1394	    }
1395
1396	    /* Store C(4,1) block. */
1397	    for (i = 0; i < 4; i++)
1398	    {
1399	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000004624 <+17956>:	mulps  %xmm6,%xmm9
   0x000000000000463f <+17983>:	mulps  %xmm6,%xmm8
   0x0000000000004657 <+18007>:	mulps  %xmm6,%xmm7
   0x000000000000466d <+18029>:	mulps  %xmm6,%xmm11

1400	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_41]), C_row[i]);
   0x0000000000004628 <+17960>:	addps  0x300(%r9),%xmm9
   0x0000000000004643 <+17987>:	addps  0x310(%r9),%xmm8
   0x000000000000465a <+18010>:	addps  0x320(%r9),%xmm7
   0x0000000000004671 <+18033>:	addps  0x330(%r9),%xmm11

1401	      _mm_store_ps(&C[i*4+C_OFFSET_41], C_row[i]);
   0x0000000000004637 <+17975>:	movaps %xmm9,0x300(%r9)
   0x000000000000464b <+17995>:	movaps %xmm8,0x310(%r9)
   0x0000000000004662 <+18018>:	movaps %xmm7,0x320(%r9)
   0x0000000000004679 <+18041>:	movaps %xmm11,0x330(%r9)

1402	    }
1403
1404	    /* Reset C(4,2) matrix accumulators */
1405	    C_row[0] = _mm_setzero_ps();
   0x0000000000004653 <+18003>:	xorps  %xmm8,%xmm8

1406	    C_row[1] = _mm_setzero_ps();
   0x000000000000466a <+18026>:	xorps  %xmm7,%xmm7

1407	    C_row[2] = _mm_setzero_ps();
   0x0000000000004630 <+17968>:	xorps  %xmm5,%xmm5

1408	    C_row[3] = _mm_setzero_ps();
   0x0000000000004633 <+17971>:	xorps  %xmm10,%xmm10

1409
1410	    if (norm[12]*norm[17] >= tolerance &&
   0x0000000000004681 <+18049>:	mulss  %xmm0,%xmm4
   0x0000000000004685 <+18053>:	comiss %xmm1,%xmm4
   0x0000000000004688 <+18056>:	jb     0x4b2b <stream_kernel+19243>

1411	        norm[13]*norm[21] >= tolerance &&
   0x000000000000468e <+18062>:	movss  0x4c(%rdx,%rsi,1),%xmm4
   0x0000000000004694 <+18068>:	mulss  0x6c(%rdx,%rsi,1),%xmm4
   0x000000000000469a <+18074>:	comiss %xmm1,%xmm4
   0x000000000000469d <+18077>:	jb     0x4b2b <stream_kernel+19243>

1412	        norm[14]*norm[25] >= tolerance &&
   0x00000000000046a3 <+18083>:	movss  0x50(%rdx,%rsi,1),%xmm4
   0x00000000000046a9 <+18089>:	mulss  0x7c(%rdx,%rsi,1),%xmm4
   0x00000000000046af <+18095>:	comiss %xmm1,%xmm4
   0x00000000000046b2 <+18098>:	jb     0x4b2b <stream_kernel+19243>

1413	        norm[15]*norm[29] >= tolerance)
   0x00000000000046b8 <+18104>:	movss  0x54(%rdx,%rsi,1),%xmm4
   0x00000000000046be <+18110>:	mulss  0x8c(%rdx,%rsi,1),%xmm4
   0x00000000000046c7 <+18119>:	comiss %xmm1,%xmm4
   0x00000000000046ca <+18122>:	jb     0x4b2b <stream_kernel+19243>

1414	    {
1415	      /* A(4,1)*B(1,2) = C(4,2). */
1416	      for (i = 0; i < 4; i++)
1417	      {
1418	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1419	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
1420	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000046d0 <+18128>:	movaps 0x40(%rcx),%xmm11
   0x00000000000046e3 <+18147>:	movaps 0xc00(%r8),%xmm7
   0x0000000000004703 <+18179>:	movaps 0xc40(%r8),%xmm8
   0x0000000000004713 <+18195>:	mulps  %xmm11,%xmm7
   0x0000000000004747 <+18247>:	movaps 0xc80(%r8),%xmm13
   0x000000000000474f <+18255>:	mulps  %xmm11,%xmm8
   0x000000000000477a <+18298>:	movaps 0xcc0(%r8),%xmm10
   0x0000000000004782 <+18306>:	mulps  %xmm11,%xmm13
   0x00000000000047b6 <+18358>:	mulps  %xmm11,%xmm10

1421
1422	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1423	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
1424	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000046d5 <+18133>:	movaps 0x50(%rcx),%xmm14
   0x00000000000046eb <+18155>:	movaps 0xc10(%r8),%xmm10
   0x000000000000470b <+18187>:	movaps 0xc50(%r8),%xmm5
   0x0000000000004717 <+18199>:	mulps  %xmm14,%xmm10
   0x0000000000004723 <+18211>:	addps  %xmm7,%xmm10
   0x0000000000004753 <+18259>:	mulps  %xmm14,%xmm5
   0x0000000000004757 <+18263>:	addps  %xmm8,%xmm5
   0x000000000000475b <+18267>:	movaps 0xc90(%r8),%xmm8
   0x0000000000004786 <+18310>:	mulps  %xmm14,%xmm8
   0x000000000000478a <+18314>:	addps  %xmm13,%xmm8
   0x00000000000047ba <+18362>:	movaps 0xcd0(%r8),%xmm11
   0x00000000000047c2 <+18370>:	mulps  %xmm14,%xmm11
   0x00000000000047ce <+18382>:	addps  %xmm10,%xmm11

1425
1426	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1427	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
1428	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000046da <+18138>:	movaps 0x60(%rcx),%xmm12
   0x00000000000046f3 <+18163>:	movaps 0xc20(%r8),%xmm13
   0x000000000000471b <+18203>:	movaps 0xca0(%r8),%xmm15
   0x000000000000472f <+18223>:	mulps  %xmm12,%xmm13
   0x0000000000004733 <+18227>:	addps  %xmm10,%xmm13
   0x0000000000004737 <+18231>:	movaps 0xc60(%r8),%xmm10
   0x0000000000004763 <+18275>:	mulps  %xmm12,%xmm10
   0x0000000000004767 <+18279>:	addps  %xmm5,%xmm10
   0x0000000000004796 <+18326>:	mulps  %xmm12,%xmm15
   0x000000000000479a <+18330>:	addps  %xmm8,%xmm15
   0x00000000000047c6 <+18374>:	movaps 0xce0(%r8),%xmm14
   0x00000000000047d2 <+18386>:	mulps  %xmm12,%xmm14
   0x00000000000047de <+18398>:	addps  %xmm11,%xmm14

1429
1430	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1431	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
1432	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000046df <+18143>:	movaps 0x70(%rcx),%xmm4
   0x00000000000046fb <+18171>:	movaps 0xc30(%r8),%xmm9
   0x0000000000004727 <+18215>:	movaps 0xc70(%r8),%xmm7
   0x000000000000473f <+18239>:	mulps  %xmm4,%xmm9
   0x0000000000004743 <+18243>:	addps  %xmm13,%xmm9
   0x000000000000476b <+18283>:	movaps 0xcf0(%r8),%xmm5
   0x0000000000004773 <+18291>:	mulps  %xmm4,%xmm7
   0x0000000000004776 <+18294>:	addps  %xmm10,%xmm7
   0x000000000000479e <+18334>:	movaps 0xcb0(%r8),%xmm8
   0x00000000000047a6 <+18342>:	mulps  %xmm4,%xmm8
   0x00000000000047aa <+18346>:	addps  %xmm15,%xmm8
   0x00000000000047ea <+18410>:	mulps  %xmm4,%xmm5
   0x00000000000047f9 <+18425>:	addps  %xmm14,%xmm5

1433	      }
1434
1435	      /* A(4,2)*B(2,2) = C(4,2). */
1436	      for (i = 0; i < 4; i++)
1437	      {
1438	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1439	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
1440	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000478e <+18318>:	movaps 0x140(%rcx),%xmm13
   0x00000000000047d6 <+18390>:	movaps 0xd00(%r8),%xmm12
   0x00000000000047e2 <+18402>:	movaps 0xd40(%r8),%xmm11
   0x00000000000047f5 <+18421>:	mulps  %xmm13,%xmm12
   0x00000000000047fd <+18429>:	mulps  %xmm13,%xmm11
   0x0000000000004815 <+18453>:	addps  %xmm9,%xmm12
   0x0000000000004835 <+18485>:	addps  %xmm7,%xmm11
   0x0000000000004849 <+18505>:	movaps 0xd80(%r8),%xmm4
   0x0000000000004851 <+18513>:	mulps  %xmm13,%xmm4
   0x0000000000004895 <+18581>:	addps  %xmm8,%xmm4
   0x00000000000048a9 <+18601>:	movaps 0xdc0(%r8),%xmm4
   0x00000000000048b1 <+18609>:	mulps  %xmm13,%xmm4
   0x00000000000048f1 <+18673>:	addps  %xmm5,%xmm4

1441
1442	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1443	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
1444	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000047ed <+18413>:	movaps 0xd10(%r8),%xmm4
   0x0000000000004801 <+18433>:	movaps 0x150(%rcx),%xmm10
   0x0000000000004811 <+18449>:	mulps  %xmm10,%xmm4
   0x0000000000004825 <+18469>:	addps  %xmm12,%xmm4
   0x0000000000004839 <+18489>:	movaps 0xd50(%r8),%xmm7
   0x0000000000004841 <+18497>:	mulps  %xmm10,%xmm7
   0x0000000000004865 <+18533>:	addps  %xmm11,%xmm7
   0x0000000000004899 <+18585>:	movaps 0xd90(%r8),%xmm8
   0x00000000000048a1 <+18593>:	mulps  %xmm10,%xmm8
   0x00000000000048a5 <+18597>:	addps  %xmm4,%xmm8
   0x00000000000048b5 <+18613>:	movaps 0xdd0(%r8),%xmm13
   0x00000000000048bd <+18621>:	mulps  %xmm10,%xmm13
   0x000000000000490c <+18700>:	addps  %xmm4,%xmm13

1445
1446	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1447	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
1448	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004809 <+18441>:	movaps 0x160(%rcx),%xmm14
   0x0000000000004829 <+18473>:	movaps 0xd20(%r8),%xmm12
   0x0000000000004831 <+18481>:	mulps  %xmm14,%xmm12
   0x0000000000004845 <+18501>:	addps  %xmm4,%xmm12
   0x0000000000004859 <+18521>:	movaps 0xd60(%r8),%xmm12
   0x0000000000004861 <+18529>:	mulps  %xmm14,%xmm12
   0x0000000000004875 <+18549>:	addps  %xmm7,%xmm12
   0x0000000000004879 <+18553>:	movaps 0xda0(%r8),%xmm7
   0x0000000000004881 <+18561>:	mulps  %xmm14,%xmm7
   0x00000000000048d5 <+18645>:	addps  %xmm8,%xmm7
   0x00000000000048f4 <+18676>:	movaps 0xde0(%r8),%xmm5
   0x00000000000048fc <+18684>:	mulps  %xmm14,%xmm5
   0x000000000000491b <+18715>:	addps  %xmm13,%xmm5

1449
1450	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1451	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
1452	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000047ae <+18350>:	movaps 0x170(%rcx),%xmm15
   0x0000000000004819 <+18457>:	movaps 0xd30(%r8),%xmm9
   0x0000000000004821 <+18465>:	mulps  %xmm15,%xmm9
   0x0000000000004855 <+18517>:	addps  %xmm12,%xmm9
   0x0000000000004869 <+18537>:	movaps 0xd70(%r8),%xmm11
   0x0000000000004871 <+18545>:	mulps  %xmm15,%xmm11
   0x0000000000004885 <+18565>:	addps  %xmm12,%xmm11
   0x0000000000004889 <+18569>:	movaps 0xdb0(%r8),%xmm12
   0x0000000000004891 <+18577>:	mulps  %xmm15,%xmm12
   0x00000000000048c1 <+18625>:	movaps 0xdf0(%r8),%xmm10
   0x00000000000048c9 <+18633>:	mulps  %xmm15,%xmm10
   0x00000000000048e1 <+18657>:	addps  %xmm7,%xmm12
   0x0000000000004927 <+18727>:	addps  %xmm5,%xmm10

1453	      }
1454
1455	      /* A(4,3)*B(3,2) = C(4,2). */
1456	      for (i = 0; i < 4; i++)
1457	      {
1458	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1459	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_32]);
1460	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000048d9 <+18649>:	movaps 0x240(%rcx),%xmm8
   0x00000000000048e5 <+18661>:	movaps 0xe40(%r8),%xmm7
   0x00000000000048ed <+18669>:	mulps  %xmm8,%xmm7
   0x0000000000004900 <+18688>:	movaps 0xe00(%r8),%xmm14
   0x0000000000004908 <+18696>:	mulps  %xmm8,%xmm14
   0x0000000000004932 <+18738>:	addps  %xmm9,%xmm14
   0x0000000000004942 <+18754>:	addps  %xmm11,%xmm7
   0x0000000000004956 <+18774>:	movaps 0xe80(%r8),%xmm14
   0x000000000000495e <+18782>:	mulps  %xmm8,%xmm14
   0x0000000000004995 <+18837>:	movaps 0xec0(%r8),%xmm11
   0x000000000000499d <+18845>:	mulps  %xmm8,%xmm11
   0x00000000000049bc <+18876>:	addps  %xmm12,%xmm14
   0x00000000000049cc <+18892>:	addps  %xmm10,%xmm11

1461
1462	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1463	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_32]);
1464	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000492b <+18731>:	movaps 0x250(%rcx),%xmm5
   0x0000000000004936 <+18742>:	movaps 0xe10(%r8),%xmm9
   0x000000000000493e <+18750>:	mulps  %xmm5,%xmm9
   0x0000000000004946 <+18758>:	movaps 0xe50(%r8),%xmm11
   0x000000000000494e <+18766>:	mulps  %xmm5,%xmm11
   0x0000000000004952 <+18770>:	addps  %xmm14,%xmm9
   0x0000000000004972 <+18802>:	addps  %xmm7,%xmm11
   0x00000000000049c0 <+18880>:	movaps 0xe90(%r8),%xmm12
   0x00000000000049c8 <+18888>:	mulps  %xmm5,%xmm12
   0x00000000000049d0 <+18896>:	movaps 0xed0(%r8),%xmm10
   0x00000000000049e0 <+18912>:	mulps  %xmm5,%xmm10
   0x00000000000049f0 <+18928>:	addps  %xmm14,%xmm12
   0x0000000000004a20 <+18976>:	addps  %xmm11,%xmm10

1465
1466	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1467	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_32]);
1468	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000048cd <+18637>:	movaps 0xe20(%r8),%xmm15
   0x0000000000004910 <+18704>:	movaps 0x260(%rcx),%xmm4
   0x0000000000004917 <+18711>:	mulps  %xmm4,%xmm15
   0x0000000000004962 <+18786>:	addps  %xmm9,%xmm15
   0x0000000000004976 <+18806>:	movaps 0xe60(%r8),%xmm7
   0x000000000000497e <+18814>:	mulps  %xmm4,%xmm7
   0x0000000000004991 <+18833>:	addps  %xmm11,%xmm7
   0x00000000000049a1 <+18849>:	movaps 0xee0(%r8),%xmm8
   0x00000000000049a9 <+18857>:	mulps  %xmm4,%xmm8
   0x00000000000049b1 <+18865>:	movaps 0xea0(%r8),%xmm7
   0x00000000000049b9 <+18873>:	mulps  %xmm4,%xmm7
   0x0000000000004a00 <+18944>:	addps  %xmm12,%xmm7
   0x0000000000004a30 <+18992>:	addps  %xmm10,%xmm8

1469
1470	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1471	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_32]);
1472	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000491f <+18719>:	movaps 0xe30(%r8),%xmm13
   0x0000000000004966 <+18790>:	movaps 0x270(%rcx),%xmm9
   0x000000000000496e <+18798>:	mulps  %xmm9,%xmm13
   0x0000000000004981 <+18817>:	addps  %xmm15,%xmm13
   0x0000000000004985 <+18821>:	movaps 0xe70(%r8),%xmm15
   0x000000000000498d <+18829>:	mulps  %xmm9,%xmm15
   0x00000000000049ad <+18861>:	addps  %xmm7,%xmm15
   0x00000000000049d8 <+18904>:	movaps 0xef0(%r8),%xmm4
   0x00000000000049ec <+18924>:	mulps  %xmm9,%xmm4
   0x00000000000049f4 <+18932>:	movaps 0xeb0(%r8),%xmm14
   0x00000000000049fc <+18940>:	mulps  %xmm9,%xmm14
   0x0000000000004a0c <+18956>:	addps  %xmm7,%xmm14
   0x0000000000004a40 <+19008>:	addps  %xmm8,%xmm4

1473	      }
1474
1475	      /* A(4,4)*B(4,2) = C(4,2). */
1476	      for (i = 0; i < 4; i++)
1477	      {
1478	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1479	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_42]);
1480	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000049e4 <+18916>:	movaps 0xf00(%r8),%xmm5
   0x0000000000004a18 <+18968>:	movaps 0xf40(%r8),%xmm9
   0x0000000000004a24 <+18980>:	movaps 0x340(%rcx),%xmm11
   0x0000000000004a2c <+18988>:	mulps  %xmm11,%xmm5
   0x0000000000004a44 <+19012>:	mulps  %xmm11,%xmm9
   0x0000000000004a50 <+19024>:	addps  %xmm13,%xmm5
   0x0000000000004a60 <+19040>:	addps  %xmm15,%xmm9
   0x0000000000004a94 <+19092>:	movaps 0xf80(%r8),%xmm9
   0x0000000000004a9c <+19100>:	mulps  %xmm11,%xmm9
   0x0000000000004ab4 <+19124>:	movaps 0xfc0(%r8),%xmm15
   0x0000000000004abc <+19132>:	mulps  %xmm11,%xmm15
   0x0000000000004acf <+19151>:	addps  %xmm14,%xmm9
   0x0000000000004adf <+19167>:	addps  %xmm4,%xmm15

1481
1482	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1483	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_42]);
1484	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004a04 <+18948>:	movaps 0x350(%rcx),%xmm12
   0x0000000000004a54 <+19028>:	movaps 0xf10(%r8),%xmm13
   0x0000000000004a5c <+19036>:	mulps  %xmm12,%xmm13
   0x0000000000004a64 <+19044>:	movaps 0xf50(%r8),%xmm15
   0x0000000000004a6c <+19052>:	mulps  %xmm12,%xmm15
   0x0000000000004a70 <+19056>:	addps  %xmm5,%xmm13
   0x0000000000004a90 <+19088>:	addps  %xmm9,%xmm15
   0x0000000000004ad3 <+19155>:	movaps 0xf90(%r8),%xmm14
   0x0000000000004adb <+19163>:	mulps  %xmm12,%xmm14
   0x0000000000004ae3 <+19171>:	movaps 0xfd0(%r8),%xmm4
   0x0000000000004aeb <+19179>:	mulps  %xmm12,%xmm4
   0x0000000000004aef <+19183>:	addps  %xmm9,%xmm14
   0x0000000000004aff <+19199>:	addps  %xmm15,%xmm4

1485
1486	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1487	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_42]);
1488	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004a10 <+18960>:	movaps 0xf20(%r8),%xmm7
   0x0000000000004a34 <+18996>:	movaps 0x360(%rcx),%xmm10
   0x0000000000004a3c <+19004>:	mulps  %xmm10,%xmm7
   0x0000000000004a74 <+19060>:	movaps 0xf60(%r8),%xmm5
   0x0000000000004a7c <+19068>:	mulps  %xmm10,%xmm5
   0x0000000000004a80 <+19072>:	addps  %xmm13,%xmm7
   0x0000000000004ab0 <+19120>:	addps  %xmm15,%xmm5
   0x0000000000004af3 <+19187>:	movaps 0xfa0(%r8),%xmm9
   0x0000000000004afb <+19195>:	mulps  %xmm10,%xmm9
   0x0000000000004b03 <+19203>:	addps  %xmm14,%xmm9
   0x0000000000004b0b <+19211>:	movaps 0xfe0(%r8),%xmm9
   0x0000000000004b13 <+19219>:	mulps  %xmm10,%xmm9
   0x0000000000004b23 <+19235>:	addps  %xmm4,%xmm9

1489
1490	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1491	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_42]);
1492	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004a48 <+19016>:	movaps 0xf30(%r8),%xmm8
   0x0000000000004a84 <+19076>:	movaps 0x370(%rcx),%xmm13
   0x0000000000004a8c <+19084>:	mulps  %xmm13,%xmm8
   0x0000000000004aa0 <+19104>:	addps  %xmm7,%xmm8
   0x0000000000004aa4 <+19108>:	movaps 0xf70(%r8),%xmm7
   0x0000000000004aac <+19116>:	mulps  %xmm13,%xmm7
   0x0000000000004ac0 <+19136>:	addps  %xmm5,%xmm7
   0x0000000000004ac3 <+19139>:	movaps 0xfb0(%r8),%xmm5
   0x0000000000004acb <+19147>:	mulps  %xmm13,%xmm5
   0x0000000000004b07 <+19207>:	addps  %xmm9,%xmm5
   0x0000000000004b17 <+19223>:	movaps 0xff0(%r8),%xmm10
   0x0000000000004b1f <+19231>:	mulps  %xmm13,%xmm10
   0x0000000000004b27 <+19239>:	addps  %xmm9,%xmm10

1493	      }
1494	    }
1495
1496	    /* Store C(4,2) block. */
1497	    for (i = 0; i < 4; i++)
1498	    {
1499	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000004b2b <+19243>:	mulps  %xmm6,%xmm8
   0x0000000000004b46 <+19270>:	mulps  %xmm6,%xmm7
   0x0000000000004b5c <+19292>:	mulps  %xmm6,%xmm5
   0x0000000000004b72 <+19314>:	mulps  %xmm6,%xmm10

1500	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_42]), C_row[i]);
   0x0000000000004b2f <+19247>:	addps  0x340(%r9),%xmm8
   0x0000000000004b49 <+19273>:	addps  0x350(%r9),%xmm7
   0x0000000000004b5f <+19295>:	addps  0x360(%r9),%xmm5
   0x0000000000004b76 <+19318>:	addps  0x370(%r9),%xmm10

1501	      _mm_store_ps(&C[i*4+C_OFFSET_42], C_row[i]);
   0x0000000000004b3a <+19258>:	movaps %xmm8,0x340(%r9)
   0x0000000000004b51 <+19281>:	movaps %xmm7,0x350(%r9)
   0x0000000000004b67 <+19303>:	movaps %xmm5,0x360(%r9)
   0x0000000000004b7e <+19326>:	movaps %xmm10,0x370(%r9)

1502	    }
1503
1504	    /* Reset C(4,3) matrix accumulators */
1505	    C_row[0] = _mm_setzero_ps();
   0x0000000000004b59 <+19289>:	xorps  %xmm7,%xmm7

1506	    C_row[1] = _mm_setzero_ps();
   0x0000000000004b6f <+19311>:	xorps  %xmm5,%xmm5

1507	    C_row[2] = _mm_setzero_ps();
   0x0000000000004b37 <+19255>:	xorps  %xmm4,%xmm4

1508	    C_row[3] = _mm_setzero_ps();
   0x0000000000004b42 <+19266>:	xorps  %xmm8,%xmm8

1509
1510	    if (norm[12]*norm[18] >= tolerance &&
   0x0000000000004b86 <+19334>:	mulss  %xmm0,%xmm3
   0x0000000000004b8a <+19338>:	comiss %xmm1,%xmm3
   0x0000000000004b8d <+19341>:	jb     0x5039 <stream_kernel+20537>

1511	        norm[13]*norm[22] >= tolerance &&
   0x0000000000004b93 <+19347>:	movss  0x4c(%rdx,%rsi,1),%xmm3
   0x0000000000004b99 <+19353>:	mulss  0x70(%rdx,%rsi,1),%xmm3
   0x0000000000004b9f <+19359>:	comiss %xmm1,%xmm3
   0x0000000000004ba2 <+19362>:	jb     0x5039 <stream_kernel+20537>

1512	        norm[14]*norm[26] >= tolerance &&
   0x0000000000004ba8 <+19368>:	movss  0x50(%rdx,%rsi,1),%xmm3
   0x0000000000004bae <+19374>:	mulss  0x80(%rdx,%rsi,1),%xmm3
   0x0000000000004bb7 <+19383>:	comiss %xmm1,%xmm3
   0x0000000000004bba <+19386>:	jb     0x5039 <stream_kernel+20537>

1513	        norm[15]*norm[30] >= tolerance)
   0x0000000000004bc0 <+19392>:	movss  0x54(%rdx,%rsi,1),%xmm3
   0x0000000000004bc6 <+19398>:	mulss  0x90(%rdx,%rsi,1),%xmm3
   0x0000000000004bcf <+19407>:	comiss %xmm1,%xmm3
   0x0000000000004bd2 <+19410>:	jb     0x5039 <stream_kernel+20537>

1514	    {
1515	      /* A(4,1)*B(1,3) = C(4,3). */
1516	      for (i = 0; i < 4; i++)
1517	      {
1518	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1519	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_13]);
1520	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004bd8 <+19416>:	movaps 0x80(%rcx),%xmm11
   0x0000000000004bef <+19439>:	movaps 0xc00(%r8),%xmm3
   0x0000000000004c0f <+19471>:	movaps 0xc40(%r8),%xmm14
   0x0000000000004c27 <+19495>:	movaps 0xc80(%r8),%xmm9
   0x0000000000004c37 <+19511>:	mulps  %xmm11,%xmm3
   0x0000000000004c67 <+19559>:	mulps  %xmm11,%xmm14
   0x0000000000004c83 <+19587>:	movaps 0xcc0(%r8),%xmm13
   0x0000000000004c9a <+19610>:	mulps  %xmm11,%xmm9
   0x0000000000004ccd <+19661>:	mulps  %xmm11,%xmm13

1521
1522	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1523	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_13]);
1524	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004be0 <+19424>:	movaps 0x90(%rcx),%xmm10
   0x0000000000004bf7 <+19447>:	movaps 0xc10(%r8),%xmm5
   0x0000000000004c17 <+19479>:	movaps 0xc50(%r8),%xmm13
   0x0000000000004c2f <+19503>:	movaps 0xc90(%r8),%xmm15
   0x0000000000004c3b <+19515>:	mulps  %xmm10,%xmm5
   0x0000000000004c3f <+19519>:	addps  %xmm3,%xmm5
   0x0000000000004c6b <+19563>:	mulps  %xmm10,%xmm13
   0x0000000000004c6f <+19567>:	addps  %xmm14,%xmm13
   0x0000000000004c9e <+19614>:	mulps  %xmm10,%xmm15
   0x0000000000004ca2 <+19618>:	addps  %xmm9,%xmm15
   0x0000000000004cd1 <+19665>:	movaps 0xcd0(%r8),%xmm11
   0x0000000000004cd9 <+19673>:	mulps  %xmm10,%xmm11
   0x0000000000004ce5 <+19685>:	addps  %xmm13,%xmm11

1525
1526	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1527	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_13]);
1528	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004be8 <+19432>:	movaps 0xa0(%rcx),%xmm7
   0x0000000000004bff <+19455>:	movaps 0xc20(%r8),%xmm4
   0x0000000000004c1f <+19487>:	movaps 0xc60(%r8),%xmm8
   0x0000000000004c4a <+19530>:	mulps  %xmm7,%xmm4
   0x0000000000004c4d <+19533>:	addps  %xmm5,%xmm4
   0x0000000000004c73 <+19571>:	movaps 0xca0(%r8),%xmm14
   0x0000000000004c7b <+19579>:	mulps  %xmm7,%xmm8
   0x0000000000004c7f <+19583>:	addps  %xmm13,%xmm8
   0x0000000000004c92 <+19602>:	movaps 0xce0(%r8),%xmm8
   0x0000000000004cae <+19630>:	mulps  %xmm7,%xmm14
   0x0000000000004cb2 <+19634>:	addps  %xmm15,%xmm14
   0x0000000000004ce9 <+19689>:	mulps  %xmm7,%xmm8
   0x0000000000004cf5 <+19701>:	addps  %xmm11,%xmm8

1529
1530	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1531	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_13]);
1532	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004c07 <+19463>:	movaps 0xc30(%r8),%xmm12
   0x0000000000004c42 <+19522>:	movaps 0xc70(%r8),%xmm3
   0x0000000000004c50 <+19536>:	movaps 0xb0(%rcx),%xmm5
   0x0000000000004c57 <+19543>:	mulps  %xmm5,%xmm12
   0x0000000000004c5b <+19547>:	addps  %xmm4,%xmm12
   0x0000000000004c5f <+19551>:	movaps 0xcb0(%r8),%xmm4
   0x0000000000004c8b <+19595>:	mulps  %xmm5,%xmm3
   0x0000000000004c8e <+19598>:	addps  %xmm8,%xmm3
   0x0000000000004cbe <+19646>:	mulps  %xmm5,%xmm4
   0x0000000000004cc1 <+19649>:	addps  %xmm14,%xmm4
   0x0000000000004cdd <+19677>:	movaps 0xcf0(%r8),%xmm10
   0x0000000000004d01 <+19713>:	mulps  %xmm5,%xmm10
   0x0000000000004d11 <+19729>:	addps  %xmm8,%xmm10

1533	      }
1534
1535	      /* A(4,2)*B(2,3) = C(4,3). */
1536	      for (i = 0; i < 4; i++)
1537	      {
1538	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1539	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_23]);
1540	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004cb6 <+19638>:	movaps 0x180(%rcx),%xmm15
   0x0000000000004ced <+19693>:	movaps 0xd00(%r8),%xmm7
   0x0000000000004d0d <+19725>:	mulps  %xmm15,%xmm7
   0x0000000000004d2d <+19757>:	addps  %xmm12,%xmm7
   0x0000000000004d35 <+19765>:	movaps 0xd40(%r8),%xmm12
   0x0000000000004d49 <+19785>:	mulps  %xmm15,%xmm12
   0x0000000000004d5d <+19805>:	addps  %xmm3,%xmm12
   0x0000000000004d65 <+19813>:	movaps 0xd80(%r8),%xmm3
   0x0000000000004d89 <+19849>:	mulps  %xmm15,%xmm3
   0x0000000000004da1 <+19873>:	movaps 0xdc0(%r8),%xmm14
   0x0000000000004dad <+19885>:	addps  %xmm4,%xmm3
   0x0000000000004dbf <+19903>:	mulps  %xmm15,%xmm14
   0x0000000000004def <+19951>:	addps  %xmm10,%xmm14

1541
1542	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1543	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_23]);
1544	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004ca6 <+19622>:	movaps 0x190(%rcx),%xmm9
   0x0000000000004cf9 <+19705>:	movaps 0xd10(%r8),%xmm11
   0x0000000000004d15 <+19733>:	mulps  %xmm9,%xmm11
   0x0000000000004d3d <+19773>:	addps  %xmm7,%xmm11
   0x0000000000004d41 <+19777>:	movaps 0xd50(%r8),%xmm7
   0x0000000000004d61 <+19809>:	mulps  %xmm9,%xmm7
   0x0000000000004d7d <+19837>:	addps  %xmm12,%xmm7
   0x0000000000004db0 <+19888>:	movaps 0xd90(%r8),%xmm4
   0x0000000000004db8 <+19896>:	mulps  %xmm9,%xmm4
   0x0000000000004dbc <+19900>:	addps  %xmm3,%xmm4
   0x0000000000004df3 <+19955>:	movaps 0xdd0(%r8),%xmm10
   0x0000000000004dfb <+19963>:	mulps  %xmm9,%xmm10
   0x0000000000004e0b <+19979>:	addps  %xmm14,%xmm10

1545
1546	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1547	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_23]);
1548	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004cc5 <+19653>:	movaps 0xd60(%r8),%xmm14
   0x0000000000004d05 <+19717>:	movaps 0xd20(%r8),%xmm5
   0x0000000000004d19 <+19737>:	movaps 0x1a0(%rcx),%xmm13
   0x0000000000004d29 <+19753>:	mulps  %xmm13,%xmm5
   0x0000000000004d31 <+19761>:	mulps  %xmm13,%xmm14
   0x0000000000004d4d <+19789>:	addps  %xmm11,%xmm5
   0x0000000000004d81 <+19841>:	movaps 0xda0(%r8),%xmm12
   0x0000000000004d8d <+19853>:	addps  %xmm7,%xmm14
   0x0000000000004d99 <+19865>:	mulps  %xmm13,%xmm12
   0x0000000000004dcb <+19915>:	addps  %xmm4,%xmm12
   0x0000000000004dcf <+19919>:	movaps 0xde0(%r8),%xmm4
   0x0000000000004dd7 <+19927>:	mulps  %xmm13,%xmm4
   0x0000000000004e27 <+20007>:	addps  %xmm10,%xmm4

1549
1550	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1551	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_23]);
1552	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004d21 <+19745>:	movaps 0x1b0(%rcx),%xmm8
   0x0000000000004d51 <+19793>:	movaps 0xd30(%r8),%xmm11
   0x0000000000004d59 <+19801>:	mulps  %xmm8,%xmm11
   0x0000000000004d6d <+19821>:	addps  %xmm5,%xmm11
   0x0000000000004d71 <+19825>:	movaps 0xd70(%r8),%xmm5
   0x0000000000004d79 <+19833>:	mulps  %xmm8,%xmm5
   0x0000000000004d91 <+19857>:	movaps 0xdb0(%r8),%xmm7
   0x0000000000004d9d <+19869>:	addps  %xmm14,%xmm5
   0x0000000000004da9 <+19881>:	mulps  %xmm8,%xmm7
   0x0000000000004de3 <+19939>:	addps  %xmm12,%xmm7
   0x0000000000004e0f <+19983>:	movaps 0xdf0(%r8),%xmm14
   0x0000000000004e17 <+19991>:	mulps  %xmm8,%xmm14
   0x0000000000004e37 <+20023>:	addps  %xmm4,%xmm14

1553	      }
1554
1555	      /* A(4,3)*B(3,3) = C(4,3). */
1556	      for (i = 0; i < 4; i++)
1557	      {
1558	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1559	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_33]);
1560	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004de7 <+19943>:	movaps 0x280(%rcx),%xmm12
   0x0000000000004dff <+19967>:	movaps 0xe00(%r8),%xmm9
   0x0000000000004e07 <+19975>:	mulps  %xmm12,%xmm9
   0x0000000000004e43 <+20035>:	addps  %xmm11,%xmm9
   0x0000000000004e73 <+20083>:	movaps 0xe40(%r8),%xmm13
   0x0000000000004e7f <+20095>:	movaps 0xe80(%r8),%xmm8
   0x0000000000004e87 <+20103>:	mulps  %xmm12,%xmm13
   0x0000000000004e8b <+20107>:	addps  %xmm5,%xmm13
   0x0000000000004eaf <+20143>:	mulps  %xmm12,%xmm8
   0x0000000000004ec7 <+20167>:	addps  %xmm7,%xmm8
   0x0000000000004eea <+20202>:	movaps 0xec0(%r8),%xmm7
   0x0000000000004ef2 <+20210>:	mulps  %xmm12,%xmm7
   0x0000000000004f0e <+20238>:	addps  %xmm14,%xmm7

1561
1562	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1563	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_33]);
1564	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004ddb <+19931>:	movaps 0xe10(%r8),%xmm13
   0x0000000000004e2b <+20011>:	movaps 0x290(%rcx),%xmm10
   0x0000000000004e33 <+20019>:	mulps  %xmm10,%xmm13
   0x0000000000004e53 <+20051>:	addps  %xmm9,%xmm13
   0x0000000000004e8f <+20111>:	movaps 0xe50(%r8),%xmm5
   0x0000000000004e97 <+20119>:	mulps  %xmm10,%xmm5
   0x0000000000004e9b <+20123>:	addps  %xmm13,%xmm5
   0x0000000000004ecb <+20171>:	movaps 0xe90(%r8),%xmm7
   0x0000000000004ed3 <+20179>:	mulps  %xmm10,%xmm7
   0x0000000000004ed7 <+20183>:	addps  %xmm8,%xmm7
   0x0000000000004ef6 <+20214>:	movaps 0xed0(%r8),%xmm12
   0x0000000000004efe <+20222>:	mulps  %xmm10,%xmm12
   0x0000000000004f2a <+20266>:	addps  %xmm7,%xmm12

1565
1566	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1567	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_33]);
1568	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004dc3 <+19907>:	movaps 0x2a0(%rcx),%xmm15
   0x0000000000004e1b <+19995>:	movaps 0xe20(%r8),%xmm8
   0x0000000000004e23 <+20003>:	mulps  %xmm15,%xmm8
   0x0000000000004e5f <+20063>:	movaps 0xea0(%r8),%xmm3
   0x0000000000004e6b <+20075>:	addps  %xmm13,%xmm8
   0x0000000000004e6f <+20079>:	mulps  %xmm15,%xmm3
   0x0000000000004e9f <+20127>:	movaps 0xe60(%r8),%xmm13
   0x0000000000004ea7 <+20135>:	mulps  %xmm15,%xmm13
   0x0000000000004eab <+20139>:	addps  %xmm5,%xmm13
   0x0000000000004ee7 <+20199>:	addps  %xmm7,%xmm3
   0x0000000000004f2e <+20270>:	movaps 0xee0(%r8),%xmm7
   0x0000000000004f36 <+20278>:	mulps  %xmm15,%xmm7
   0x0000000000004f5a <+20314>:	addps  %xmm12,%xmm7

1569
1570	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1571	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_33]);
1572	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004e3b <+20027>:	movaps 0xe30(%r8),%xmm4
   0x0000000000004e47 <+20039>:	movaps 0x2b0(%rcx),%xmm11
   0x0000000000004e4f <+20047>:	mulps  %xmm11,%xmm4
   0x0000000000004e57 <+20055>:	movaps 0xe70(%r8),%xmm9
   0x0000000000004e67 <+20071>:	mulps  %xmm11,%xmm9
   0x0000000000004e7b <+20091>:	addps  %xmm8,%xmm4
   0x0000000000004ebb <+20155>:	addps  %xmm13,%xmm9
   0x0000000000004edb <+20187>:	movaps 0xeb0(%r8),%xmm8
   0x0000000000004ee3 <+20195>:	mulps  %xmm11,%xmm8
   0x0000000000004f02 <+20226>:	addps  %xmm3,%xmm8
   0x0000000000004f12 <+20242>:	movaps 0xef0(%r8),%xmm14
   0x0000000000004f1a <+20250>:	mulps  %xmm11,%xmm14
   0x0000000000004f6a <+20330>:	addps  %xmm7,%xmm14

1573	      }
1574
1575	      /* A(4,4)*B(4,3) = C(4,3). */
1576	      for (i = 0; i < 4; i++)
1577	      {
1578	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1579	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_43]);
1580	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004ebf <+20159>:	movaps 0x380(%rcx),%xmm13
   0x0000000000004f1e <+20254>:	movaps 0xf00(%r8),%xmm11
   0x0000000000004f26 <+20262>:	mulps  %xmm13,%xmm11
   0x0000000000004f42 <+20290>:	addps  %xmm4,%xmm11
   0x0000000000004f46 <+20294>:	movaps 0xf40(%r8),%xmm4
   0x0000000000004f56 <+20310>:	mulps  %xmm13,%xmm4
   0x0000000000004f8a <+20362>:	addps  %xmm9,%xmm4
   0x0000000000004f92 <+20370>:	movaps 0xf80(%r8),%xmm9
   0x0000000000004f9e <+20382>:	mulps  %xmm13,%xmm9
   0x0000000000004fb1 <+20401>:	addps  %xmm8,%xmm9
   0x0000000000004fe3 <+20451>:	movaps 0xfc0(%r8),%xmm9
   0x0000000000004feb <+20459>:	mulps  %xmm13,%xmm9
   0x0000000000004ffe <+20478>:	addps  %xmm14,%xmm9

1581
1582	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1583	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_43]);
1584	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004f3a <+20282>:	movaps 0xf10(%r8),%xmm15
   0x0000000000004f5e <+20318>:	movaps 0x390(%rcx),%xmm12
   0x0000000000004f66 <+20326>:	mulps  %xmm12,%xmm15
   0x0000000000004f7a <+20346>:	addps  %xmm11,%xmm15
   0x0000000000004fa5 <+20389>:	movaps 0xf50(%r8),%xmm5
   0x0000000000004fad <+20397>:	mulps  %xmm12,%xmm5
   0x0000000000004fb5 <+20405>:	movaps 0xf90(%r8),%xmm8
   0x0000000000004fbd <+20413>:	mulps  %xmm12,%xmm8
   0x0000000000004fc1 <+20417>:	addps  %xmm4,%xmm5
   0x0000000000004fdf <+20447>:	addps  %xmm9,%xmm8
   0x0000000000005002 <+20482>:	movaps 0xfd0(%r8),%xmm14
   0x000000000000500a <+20490>:	mulps  %xmm12,%xmm14
   0x000000000000502d <+20525>:	addps  %xmm9,%xmm14

1585
1586	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1587	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_43]);
1588	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004eb3 <+20147>:	movaps 0xf20(%r8),%xmm5
   0x0000000000004f4e <+20302>:	movaps 0xf60(%r8),%xmm3
   0x0000000000004f7e <+20350>:	movaps 0x3a0(%rcx),%xmm11
   0x0000000000004f86 <+20358>:	mulps  %xmm11,%xmm5
   0x0000000000004f8e <+20366>:	mulps  %xmm11,%xmm3
   0x0000000000004f9a <+20378>:	addps  %xmm15,%xmm5
   0x0000000000004fd0 <+20432>:	addps  %xmm5,%xmm3
   0x0000000000004ff2 <+20466>:	movaps 0xfa0(%r8),%xmm3
   0x0000000000004ffa <+20474>:	mulps  %xmm11,%xmm3
   0x000000000000500e <+20494>:	movaps 0xfe0(%r8),%xmm12
   0x0000000000005016 <+20502>:	addps  %xmm8,%xmm3
   0x0000000000005022 <+20514>:	mulps  %xmm11,%xmm12
   0x0000000000005031 <+20529>:	addps  %xmm14,%xmm12

1589
1590	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1591	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_43]);
1592	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000004f06 <+20230>:	movaps 0x3b0(%rcx),%xmm10
   0x0000000000004f6e <+20334>:	movaps 0xf30(%r8),%xmm7
   0x0000000000004f76 <+20342>:	mulps  %xmm10,%xmm7
   0x0000000000004fa2 <+20386>:	addps  %xmm5,%xmm7
   0x0000000000004fc4 <+20420>:	movaps 0xfb0(%r8),%xmm4
   0x0000000000004fcc <+20428>:	mulps  %xmm10,%xmm4
   0x0000000000004fd3 <+20435>:	movaps 0xf70(%r8),%xmm5
   0x0000000000004fdb <+20443>:	mulps  %xmm10,%xmm5
   0x0000000000004fef <+20463>:	addps  %xmm3,%xmm5
   0x000000000000501a <+20506>:	movaps 0xff0(%r8),%xmm8
   0x0000000000005026 <+20518>:	addps  %xmm3,%xmm4
   0x0000000000005029 <+20521>:	mulps  %xmm10,%xmm8
   0x0000000000005035 <+20533>:	addps  %xmm12,%xmm8

1593	      }
1594	    }
1595
1596	    /* Store C(4,3) block. */
1597	    for (i = 0; i < 4; i++)
1598	    {
1599	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000005039 <+20537>:	mulps  %xmm6,%xmm7
   0x0000000000005056 <+20566>:	mulps  %xmm6,%xmm5
   0x0000000000005069 <+20585>:	mulps  %xmm6,%xmm4
   0x000000000000507f <+20607>:	mulps  %xmm6,%xmm8

1600	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_43]), C_row[i]);
   0x000000000000503c <+20540>:	addps  0x380(%r9),%xmm7
   0x0000000000005059 <+20569>:	addps  0x390(%r9),%xmm5
   0x000000000000506c <+20588>:	addps  0x3a0(%r9),%xmm4
   0x0000000000005083 <+20611>:	addps  0x3b0(%r9),%xmm8

1601	      _mm_store_ps(&C[i*4+C_OFFSET_43], C_row[i]);
   0x000000000000504b <+20555>:	movaps %xmm7,0x380(%r9)
   0x0000000000005061 <+20577>:	movaps %xmm5,0x390(%r9)
   0x0000000000005074 <+20596>:	movaps %xmm4,0x3a0(%r9)
   0x000000000000508b <+20619>:	movaps %xmm8,0x3b0(%r9)

1602	    }
1603
1604	    /* Reset C(4,4) matrix accumulators */
1605	    C_row[0] = _mm_setzero_ps();
   0x0000000000005044 <+20548>:	xorps  %xmm11,%xmm11

1606	    C_row[1] = _mm_setzero_ps();
   0x000000000000507c <+20604>:	xorps  %xmm4,%xmm4

1607	    C_row[2] = _mm_setzero_ps();
   0x0000000000005048 <+20552>:	xorps  %xmm3,%xmm3

1608	    C_row[3] = _mm_setzero_ps();
   0x0000000000005053 <+20563>:	xorps  %xmm7,%xmm7

1609
1610	    if (norm[12]*norm[19] >= tolerance &&
   0x0000000000005093 <+20627>:	mulss  %xmm0,%xmm2
   0x0000000000005097 <+20631>:	comiss %xmm1,%xmm2
   0x000000000000509a <+20634>:	jb     0x5538 <stream_kernel+21816>

1611	        norm[13]*norm[23] >= tolerance &&
   0x00000000000050a0 <+20640>:	movss  0x4c(%rdx,%rsi,1),%xmm0
   0x00000000000050a6 <+20646>:	mulss  0x74(%rdx,%rsi,1),%xmm0
   0x00000000000050ac <+20652>:	comiss %xmm1,%xmm0
   0x00000000000050af <+20655>:	jb     0x5538 <stream_kernel+21816>

1612	        norm[14]*norm[27] >= tolerance &&
   0x00000000000050b5 <+20661>:	movss  0x50(%rdx,%rsi,1),%xmm0
   0x00000000000050bb <+20667>:	mulss  0x84(%rdx,%rsi,1),%xmm0
   0x00000000000050c4 <+20676>:	comiss %xmm1,%xmm0
   0x00000000000050c7 <+20679>:	jb     0x5538 <stream_kernel+21816>

1613	        norm[15]*norm[31] >= tolerance)
   0x00000000000050cd <+20685>:	movss  0x54(%rdx,%rsi,1),%xmm0
   0x00000000000050d3 <+20691>:	mulss  0x94(%rdx,%rsi,1),%xmm0
   0x00000000000050dc <+20700>:	comiss %xmm1,%xmm0
   0x00000000000050df <+20703>:	jb     0x5538 <stream_kernel+21816>

1614	    {
1615	      /* A(4,1)*B(1,4) = C(4,4). */
1616	      for (i = 0; i < 4; i++)
1617	      {
1618	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_41]);
1619	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_14]);
1620	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000050e5 <+20709>:	movaps 0xc0(%rcx),%xmm12
   0x00000000000050fc <+20732>:	movaps 0xc00(%r8),%xmm5
   0x000000000000511c <+20764>:	movaps 0xc40(%r8),%xmm10
   0x000000000000513c <+20796>:	movaps 0xc80(%r8),%xmm0
   0x0000000000005154 <+20820>:	mulps  %xmm12,%xmm5
   0x000000000000515f <+20831>:	movaps 0xcc0(%r8),%xmm5
   0x0000000000005186 <+20870>:	mulps  %xmm12,%xmm10
   0x00000000000051b6 <+20918>:	mulps  %xmm12,%xmm0
   0x00000000000051ea <+20970>:	mulps  %xmm12,%xmm5

1621
1622	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_41]);
1623	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_14]);
1624	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000050ed <+20717>:	movaps 0xd0(%rcx),%xmm11
   0x0000000000005104 <+20740>:	movaps 0xc10(%r8),%xmm3
   0x0000000000005124 <+20772>:	movaps 0xc50(%r8),%xmm9
   0x0000000000005144 <+20804>:	movaps 0xc90(%r8),%xmm14
   0x0000000000005158 <+20824>:	mulps  %xmm11,%xmm3
   0x000000000000515c <+20828>:	addps  %xmm5,%xmm3
   0x000000000000518a <+20874>:	mulps  %xmm11,%xmm9
   0x000000000000518e <+20878>:	addps  %xmm10,%xmm9
   0x00000000000051ba <+20922>:	mulps  %xmm11,%xmm14
   0x00000000000051be <+20926>:	addps  %xmm0,%xmm14
   0x00000000000051ee <+20974>:	movaps 0xcd0(%r8),%xmm12
   0x00000000000051f6 <+20982>:	mulps  %xmm11,%xmm12
   0x0000000000005202 <+20994>:	addps  %xmm5,%xmm12

1625
1626	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_41]);
1627	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_14]);
1628	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000050f5 <+20725>:	movaps 0xe0(%rcx),%xmm4
   0x000000000000510c <+20748>:	movaps 0xc20(%r8),%xmm8
   0x000000000000512c <+20780>:	movaps 0xc60(%r8),%xmm2
   0x000000000000514c <+20812>:	movaps 0xca0(%r8),%xmm15
   0x0000000000005167 <+20839>:	mulps  %xmm4,%xmm8
   0x000000000000516b <+20843>:	addps  %xmm3,%xmm8
   0x000000000000519a <+20890>:	mulps  %xmm4,%xmm2
   0x000000000000519d <+20893>:	addps  %xmm9,%xmm2
   0x00000000000051ca <+20938>:	mulps  %xmm4,%xmm15
   0x00000000000051ce <+20942>:	addps  %xmm14,%xmm15
   0x00000000000051fa <+20986>:	movaps 0xce0(%r8),%xmm11
   0x0000000000005206 <+20998>:	mulps  %xmm4,%xmm11
   0x0000000000005212 <+21010>:	addps  %xmm12,%xmm11

1629
1630	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_41]);
1631	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_14]);
1632	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005114 <+20756>:	movaps 0xc30(%r8),%xmm13
   0x0000000000005134 <+20788>:	movaps 0xc70(%r8),%xmm7
   0x000000000000516f <+20847>:	movaps 0xf0(%rcx),%xmm3
   0x0000000000005176 <+20854>:	mulps  %xmm3,%xmm13
   0x000000000000517a <+20858>:	addps  %xmm8,%xmm13
   0x00000000000051a9 <+20905>:	mulps  %xmm3,%xmm7
   0x00000000000051ac <+20908>:	addps  %xmm2,%xmm7
   0x00000000000051c2 <+20930>:	movaps 0xcf0(%r8),%xmm0
   0x00000000000051d2 <+20946>:	movaps 0xcb0(%r8),%xmm14
   0x00000000000051da <+20954>:	mulps  %xmm3,%xmm14
   0x00000000000051de <+20958>:	addps  %xmm15,%xmm14
   0x000000000000521e <+21022>:	mulps  %xmm3,%xmm0
   0x000000000000522c <+21036>:	addps  %xmm11,%xmm0

1633	      }
1634
1635	      /* A(4,2)*B(2,4) = C(4,4). */
1636	      for (i = 0; i < 4; i++)
1637	      {
1638	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_42]);
1639	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_24]);
1640	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000051af <+20911>:	movaps 0x1c0(%rcx),%xmm2
   0x00000000000051e2 <+20962>:	movaps 0xd80(%r8),%xmm15
   0x000000000000520a <+21002>:	movaps 0xd00(%r8),%xmm4
   0x0000000000005216 <+21014>:	movaps 0xd40(%r8),%xmm12
   0x0000000000005229 <+21033>:	mulps  %xmm2,%xmm4
   0x000000000000523c <+21052>:	addps  %xmm13,%xmm4
   0x000000000000525e <+21086>:	mulps  %xmm2,%xmm12
   0x000000000000526a <+21098>:	movaps 0xdc0(%r8),%xmm11
   0x000000000000527a <+21114>:	mulps  %xmm2,%xmm15
   0x0000000000005286 <+21126>:	addps  %xmm7,%xmm12
   0x0000000000005296 <+21142>:	addps  %xmm14,%xmm15
   0x00000000000052e1 <+21217>:	mulps  %xmm2,%xmm11
   0x0000000000005330 <+21296>:	addps  %xmm0,%xmm11

1641
1642	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_42]);
1643	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_24]);
1644	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005192 <+20882>:	movaps 0x1d0(%rcx),%xmm10
   0x0000000000005230 <+21040>:	movaps 0xd10(%r8),%xmm11
   0x0000000000005238 <+21048>:	mulps  %xmm10,%xmm11
   0x000000000000524b <+21067>:	addps  %xmm4,%xmm11
   0x000000000000528a <+21130>:	movaps 0xd50(%r8),%xmm7
   0x0000000000005292 <+21138>:	mulps  %xmm10,%xmm7
   0x000000000000529a <+21146>:	movaps 0xd90(%r8),%xmm14
   0x00000000000052a6 <+21158>:	addps  %xmm12,%xmm7
   0x00000000000052aa <+21162>:	mulps  %xmm10,%xmm14
   0x00000000000052d5 <+21205>:	addps  %xmm15,%xmm14
   0x00000000000052e9 <+21225>:	movaps 0xdd0(%r8),%xmm2
   0x00000000000052f1 <+21233>:	mulps  %xmm10,%xmm2
   0x0000000000005340 <+21312>:	addps  %xmm11,%xmm2

1645
1646	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_42]);
1647	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_24]);
1648	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005221 <+21025>:	movaps 0xd20(%r8),%xmm3
   0x0000000000005244 <+21060>:	movaps 0x1e0(%rcx),%xmm5
   0x000000000000524f <+21071>:	mulps  %xmm5,%xmm3
   0x0000000000005252 <+21074>:	movaps 0xd60(%r8),%xmm13
   0x000000000000525a <+21082>:	addps  %xmm11,%xmm3
   0x0000000000005272 <+21106>:	mulps  %xmm5,%xmm13
   0x00000000000052b6 <+21174>:	addps  %xmm7,%xmm13
   0x00000000000052be <+21182>:	movaps 0xda0(%r8),%xmm7
   0x00000000000052ca <+21194>:	mulps  %xmm5,%xmm7
   0x00000000000052e5 <+21221>:	addps  %xmm14,%xmm7
   0x00000000000052f5 <+21237>:	movaps 0xde0(%r8),%xmm10
   0x00000000000052fd <+21245>:	mulps  %xmm5,%xmm10
   0x0000000000005350 <+21328>:	addps  %xmm2,%xmm10

1649
1650	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_42]);
1651	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_24]);
1652	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000517e <+20862>:	movaps 0x1f0(%rcx),%xmm8
   0x00000000000051a1 <+20897>:	movaps 0xd30(%r8),%xmm9
   0x0000000000005240 <+21056>:	mulps  %xmm8,%xmm9
   0x0000000000005262 <+21090>:	movaps 0xdb0(%r8),%xmm4
   0x0000000000005276 <+21110>:	addps  %xmm3,%xmm9
   0x000000000000527e <+21118>:	movaps 0xd70(%r8),%xmm3
   0x00000000000052a2 <+21154>:	mulps  %xmm8,%xmm3
   0x00000000000052ba <+21178>:	mulps  %xmm8,%xmm4
   0x00000000000052c6 <+21190>:	addps  %xmm13,%xmm3
   0x0000000000005315 <+21269>:	addps  %xmm7,%xmm4
   0x0000000000005318 <+21272>:	movaps 0xdf0(%r8),%xmm7
   0x0000000000005320 <+21280>:	mulps  %xmm8,%xmm7
   0x0000000000005360 <+21344>:	addps  %xmm10,%xmm7

1653	      }
1654
1655	      /* A(4,3)*B(3,4) = C(4,4). */
1656	      for (i = 0; i < 4; i++)
1657	      {
1658	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_43]);
1659	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_34]);
1660	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000052d9 <+21209>:	movaps 0x2c0(%rcx),%xmm15
   0x0000000000005301 <+21249>:	movaps 0xe00(%r8),%xmm5
   0x0000000000005311 <+21265>:	mulps  %xmm15,%xmm5
   0x0000000000005364 <+21348>:	movaps 0xe80(%r8),%xmm10
   0x000000000000536c <+21356>:	mulps  %xmm15,%xmm10
   0x0000000000005370 <+21360>:	addps  %xmm9,%xmm5
   0x0000000000005384 <+21380>:	movaps 0xe40(%r8),%xmm5
   0x000000000000538c <+21388>:	mulps  %xmm15,%xmm5
   0x00000000000053ac <+21420>:	addps  %xmm3,%xmm5
   0x00000000000053bb <+21435>:	addps  %xmm4,%xmm10
   0x00000000000053de <+21470>:	movaps 0xec0(%r8),%xmm3
   0x00000000000053e6 <+21478>:	mulps  %xmm15,%xmm3
   0x0000000000005433 <+21555>:	addps  %xmm7,%xmm3

1661
1662	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_43]);
1663	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_34]);
1664	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000052cd <+21197>:	movaps 0x2d0(%rcx),%xmm13
   0x0000000000005324 <+21284>:	movaps 0xe10(%r8),%xmm8
   0x000000000000532c <+21292>:	mulps  %xmm13,%xmm8
   0x0000000000005380 <+21376>:	addps  %xmm5,%xmm8
   0x00000000000053af <+21423>:	movaps 0xe50(%r8),%xmm3
   0x00000000000053b7 <+21431>:	mulps  %xmm13,%xmm3
   0x00000000000053bf <+21439>:	movaps 0xe90(%r8),%xmm4
   0x00000000000053c7 <+21447>:	mulps  %xmm13,%xmm4
   0x00000000000053cb <+21451>:	addps  %xmm5,%xmm3
   0x00000000000053f6 <+21494>:	addps  %xmm10,%xmm4
   0x0000000000005436 <+21558>:	movaps 0xed0(%r8),%xmm7
   0x000000000000543e <+21566>:	mulps  %xmm13,%xmm7
   0x000000000000545a <+21594>:	addps  %xmm3,%xmm7

1665
1666	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_43]);
1667	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_34]);
1668	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005348 <+21320>:	movaps 0x2e0(%rcx),%xmm11
   0x0000000000005374 <+21364>:	movaps 0xe20(%r8),%xmm9
   0x000000000000537c <+21372>:	mulps  %xmm11,%xmm9
   0x0000000000005390 <+21392>:	addps  %xmm8,%xmm9
   0x0000000000005394 <+21396>:	movaps 0xe60(%r8),%xmm8
   0x000000000000539c <+21404>:	mulps  %xmm11,%xmm8
   0x00000000000053ce <+21454>:	movaps 0xea0(%r8),%xmm5
   0x00000000000053d6 <+21462>:	mulps  %xmm11,%xmm5
   0x00000000000053da <+21466>:	addps  %xmm3,%xmm8
   0x000000000000541a <+21530>:	addps  %xmm4,%xmm5
   0x0000000000005442 <+21570>:	movaps 0xee0(%r8),%xmm13
   0x000000000000544a <+21578>:	mulps  %xmm11,%xmm13
   0x0000000000005465 <+21605>:	addps  %xmm7,%xmm13

1669
1670	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_43]);
1671	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_34]);
1672	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000052ae <+21166>:	movaps 0x2f0(%rcx),%xmm12
   0x0000000000005309 <+21257>:	movaps 0xe30(%r8),%xmm14
   0x0000000000005334 <+21300>:	mulps  %xmm12,%xmm14
   0x0000000000005338 <+21304>:	movaps 0xe70(%r8),%xmm0
   0x0000000000005344 <+21316>:	mulps  %xmm12,%xmm0
   0x0000000000005354 <+21332>:	movaps 0xeb0(%r8),%xmm2
   0x000000000000535c <+21340>:	mulps  %xmm12,%xmm2
   0x00000000000053a0 <+21408>:	addps  %xmm9,%xmm14
   0x00000000000053ea <+21482>:	addps  %xmm8,%xmm0
   0x00000000000053fa <+21498>:	movaps 0xef0(%r8),%xmm10
   0x000000000000540a <+21514>:	mulps  %xmm12,%xmm10
   0x0000000000005429 <+21545>:	addps  %xmm5,%xmm2
   0x0000000000005474 <+21620>:	addps  %xmm13,%xmm10

1673	      }
1674
1675	      /* A(4,4)*B(4,4) = C(4,4). */
1676	      for (i = 0; i < 4; i++)
1677	      {
1678	        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_44]);
1679	        B_row = _mm_load_ps(&B[0*4+B_OFFSET_44]);
1680	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000053a4 <+21412>:	movaps 0x3c0(%rcx),%xmm9
   0x000000000000541d <+21533>:	movaps 0xf40(%r8),%xmm4
   0x0000000000005425 <+21541>:	mulps  %xmm9,%xmm4
   0x000000000000544e <+21582>:	movaps 0xf00(%r8),%xmm11
   0x0000000000005456 <+21590>:	mulps  %xmm9,%xmm11
   0x0000000000005483 <+21635>:	addps  %xmm14,%xmm11
   0x0000000000005487 <+21639>:	movaps 0xfc0(%r8),%xmm14
   0x0000000000005497 <+21655>:	mulps  %xmm9,%xmm14
   0x00000000000054ab <+21675>:	movaps 0xf80(%r8),%xmm12
   0x00000000000054b3 <+21683>:	addps  %xmm0,%xmm4
   0x00000000000054c6 <+21702>:	mulps  %xmm9,%xmm12
   0x00000000000054e7 <+21735>:	addps  %xmm2,%xmm12
   0x0000000000005528 <+21800>:	addps  %xmm10,%xmm14

1681
1682	        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_44]);
1683	        B_row = _mm_load_ps(&B[1*4+B_OFFSET_44]);
1684	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x00000000000053ee <+21486>:	movaps 0x3d0(%rcx),%xmm8
   0x000000000000540e <+21518>:	movaps 0xf10(%r8),%xmm12
   0x0000000000005416 <+21526>:	mulps  %xmm8,%xmm12
   0x0000000000005493 <+21651>:	addps  %xmm11,%xmm12
   0x00000000000054b6 <+21686>:	movaps 0xf50(%r8),%xmm0
   0x00000000000054be <+21694>:	mulps  %xmm8,%xmm0
   0x00000000000054ca <+21706>:	movaps 0xfd0(%r8),%xmm9
   0x00000000000054d2 <+21714>:	addps  %xmm4,%xmm0
   0x00000000000054d5 <+21717>:	mulps  %xmm8,%xmm9
   0x00000000000054eb <+21739>:	movaps 0xf90(%r8),%xmm2
   0x00000000000054f3 <+21747>:	mulps  %xmm8,%xmm2
   0x000000000000550d <+21773>:	addps  %xmm12,%xmm2
   0x000000000000552c <+21804>:	addps  %xmm14,%xmm9

1685
1686	        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_44]);
1687	        B_row = _mm_load_ps(&B[2*4+B_OFFSET_44]);
1688	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x0000000000005402 <+21506>:	movaps 0xf20(%r8),%xmm15
   0x000000000000545d <+21597>:	movaps 0xf60(%r8),%xmm3
   0x0000000000005469 <+21609>:	movaps 0x3e0(%rcx),%xmm7
   0x0000000000005470 <+21616>:	mulps  %xmm7,%xmm15
   0x0000000000005478 <+21624>:	mulps  %xmm7,%xmm3
   0x000000000000547b <+21627>:	movaps 0xfa0(%r8),%xmm13
   0x000000000000548f <+21647>:	mulps  %xmm7,%xmm13
   0x00000000000054a3 <+21667>:	addps  %xmm12,%xmm15
   0x00000000000054e1 <+21729>:	addps  %xmm0,%xmm3
   0x00000000000054f7 <+21751>:	movaps 0xfe0(%r8),%xmm8
   0x0000000000005511 <+21777>:	mulps  %xmm7,%xmm8
   0x000000000000551d <+21789>:	addps  %xmm2,%xmm13
   0x0000000000005530 <+21808>:	addps  %xmm9,%xmm8

1689
1690	        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_44]);
1691	        B_row = _mm_load_ps(&B[3*4+B_OFFSET_44]);
1692	        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
   0x000000000000542c <+21548>:	movaps 0x3f0(%rcx),%xmm5
   0x000000000000549b <+21659>:	movaps 0xf30(%r8),%xmm11
   0x00000000000054a7 <+21671>:	mulps  %xmm5,%xmm11
   0x00000000000054c2 <+21698>:	addps  %xmm15,%xmm11
   0x00000000000054d9 <+21721>:	movaps 0xf70(%r8),%xmm4
   0x00000000000054e4 <+21732>:	mulps  %xmm5,%xmm4
   0x00000000000054ff <+21759>:	addps  %xmm3,%xmm4
   0x0000000000005502 <+21762>:	movaps 0xfb0(%r8),%xmm3
   0x000000000000550a <+21770>:	mulps  %xmm5,%xmm3
   0x0000000000005515 <+21781>:	movaps 0xff0(%r8),%xmm7
   0x0000000000005521 <+21793>:	mulps  %xmm5,%xmm7
   0x0000000000005524 <+21796>:	addps  %xmm13,%xmm3
   0x0000000000005534 <+21812>:	addps  %xmm8,%xmm7

1693	      }
1694	    }
1695
1696	    /* Store C(4,4) block. */
1697	    for (i = 0; i < 4; i++)
1698	    {
1699	      C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
   0x0000000000005538 <+21816>:	mulps  %xmm6,%xmm11
   0x0000000000005554 <+21844>:	mulps  %xmm6,%xmm4
   0x0000000000005567 <+21863>:	mulps  %xmm6,%xmm3
   0x000000000000557a <+21882>:	mulps  %xmm6,%xmm7

1700	      C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_44]), C_row[i]);
   0x000000000000553c <+21820>:	addps  0x3c0(%r9),%xmm11
   0x0000000000005557 <+21847>:	addps  0x3d0(%r9),%xmm4
   0x000000000000556a <+21866>:	addps  0x3e0(%r9),%xmm3
   0x000000000000557d <+21885>:	addps  0x3f0(%r9),%xmm7

1701	      _mm_store_ps(&C[i*4+C_OFFSET_44], C_row[i]);
   0x0000000000005548 <+21832>:	movaps %xmm11,0x3c0(%r9)
   0x000000000000555f <+21855>:	movaps %xmm4,0x3d0(%r9)
   0x0000000000005572 <+21874>:	movaps %xmm3,0x3e0(%r9)
   0x0000000000005585 <+21893>:	movaps %xmm7,0x3f0(%r9)

1702	    }
1703	  }
1704	}
   0x0000000000005593 <+21907>:	retq
   0x0000000000005594 <+21908>:	nopl   0x0(%rax,%rax,1)
   0x0000000000005599 <+21913>:	nopl   0x0(%rax)

End of assembler dump.
