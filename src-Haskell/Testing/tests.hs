import MatrixPair ; import SpAMM ; import System.Directory (removeFile)

main = do testReadWrite ; testAdd ; testMultiply

testReadWrite = do pair1 <- readPairFromMatrixMarket "Matrices/rect4arr.txt"
                   writePairToMatrixMarket pair1 "coordinate" "Matrices/tempfile.txt"
                   pair2 <- readPairFromMatrixMarket "Matrices/tempfile.txt"
                   testReport "readWrite" (pair1 == pair2)
                   removeFile "Matrices/tempfile.txt"

testAdd = do addZeros ; addZeroToValue ; addZeroToRect;
             addValueToValue ; addRectToRect

addZeros = do zeroPair <- readPairFromMatrixMarket "Matrices/zerorectcoor.txt"
              testReport "addZeros" (zeroPair `pairAdd` zeroPair == zeroPair)

addZeroToValue = do valuePair <- readPairFromMatrixMarket "Matrices/value1coor.txt"
                    zeroPair <- readPairFromMatrixMarket "Matrices/zerovaluecoor.txt"
                    testReport "addValuetoZero" (valuePair `pairAdd` zeroPair == valuePair)
                    testReport "addZerotoValue" (zeroPair `pairAdd` valuePair == valuePair)

addZeroToRect = do rectPair <- readPairFromMatrixMarket "Matrices/rect1coor.txt"
                   zeroPair <- readPairFromMatrixMarket "Matrices/zerorectcoor.txt"
                   testReport "addRectToZero" (rectPair `pairAdd` zeroPair == rectPair)
                   testReport "addZeroToRect" (zeroPair `pairAdd` rectPair == rectPair)

addValueToValue = do valuePair1 <- readPairFromMatrixMarket "Matrices/value1coor.txt"
                     valuePair2 <- readPairFromMatrixMarket "Matrices/value2coor.txt"
                     valuePair3 <- readPairFromMatrixMarket "Matrices/value3coor.txt"
                     testReport "addValueToValue" (valuePair1 `pairAdd` valuePair2 == valuePair3)

addRectToRect = do rectPair1 <- readPairFromMatrixMarket "Matrices/rect1coor.txt"
                   rectPair2 <- readPairFromMatrixMarket "Matrices/rect2coor.txt"
                   rectPair3 <- readPairFromMatrixMarket "Matrices/rect3coor.txt"
                   testReport "addRectToRect" (rectPair1 `pairAdd` rectPair2 == rectPair3)

testMultiply = do multZeros ; multZeroByValue ; multZeroByRect ; multValueByValue ;
                  multValueByRect ; multRectByRect ; multRectByRectArray

multZeros = do zeroPair <- readPairFromMatrixMarket "Matrices/zerorectcoor.txt"
               zeroSquare <- readPairFromMatrixMarket "Matrices/zerosquarecoor.txt"
               testReport "multZeros" (zeroPair `pairMult` (pairTranspose zeroPair) == zeroSquare)

multZeroByValue = do zeroRow <- readPairFromMatrixMarket "Matrices/zerorowcoor.txt"
                     valuePair <- readPairFromMatrixMarket "Matrices/value1coor.txt"
                     testReport "multValueByZero" (valuePair `pairMult` zeroRow == zeroRow)
                     testReport "multZeroByValue" ((pairTranspose zeroRow) `pairMult`
                                                  (pairTranspose valuePair) == pairTranspose zeroRow)

multZeroByRect = do zeroPair <- readPairFromMatrixMarket "Matrices/zerorectcoor.txt"
                    rectPair <- readPairFromMatrixMarket "Matrices/rect1coor.txt"
                    zeroSquare <- readPairFromMatrixMarket "Matrices/zerosquarecoor.txt"
                    testReport "multZeroByRect" (zeroPair `pairMult`
                                                (pairTranspose rectPair) == zeroSquare)
                    testReport "multZeroByRect" (rectPair `pairMult` (pairTranspose zeroPair)
                                                 == zeroSquare)

multValueByValue = do valuePair1 <- readPairFromMatrixMarket "Matrices/value2coor.txt"
                      valuePair2 <- readPairFromMatrixMarket "Matrices/value3coor.txt"
                      valuePair3 <- readPairFromMatrixMarket "Matrices/value4coor.txt"
                      testReport "multValueByValue" (valuePair1 `pairMult` valuePair2 == valuePair3)

multValueByRect = do valuePair <- readPairFromMatrixMarket "Matrices/value2coor.txt"
                     rowPair <- readPairFromMatrixMarket "Matrices/rowcoor.txt"
                     valuetimesrowPair <- readPairFromMatrixMarket "Matrices/valuetimesrowcoor.txt"
                     testReport "multValueByRect" (valuePair `pairMult` rowPair == valuetimesrowPair)
                     testReport "multRectByValue" ((pairTranspose rowPair) `pairMult`
                                                   (pairTranspose valuePair) ==
                                                    pairTranspose valuetimesrowPair)

multRectByRect = do rectPair1 <- readPairFromMatrixMarket "Matrices/rect1coor.txt"
                    rectPair2 <- readPairFromMatrixMarket "Matrices/rect2coor.txt"
                    rectPair4 <- readPairFromMatrixMarket "Matrices/rect4coor.txt"
                    testReport "multRectByRect" (rectPair1 `pairMult`
                                                (pairTranspose rectPair2) == rectPair4)

multRectByRectArray = do rectPair1 <- readPairFromMatrixMarket "Matrices/rect1arr.txt"
                         rectPair2 <- readPairFromMatrixMarket "Matrices/rect2arr.txt"
                         rectPair4 <- readPairFromMatrixMarket "Matrices/rect4arr.txt"
                         testReport "multRectByRectArray" (rectPair1 `pairMult`
                                                          (pairTranspose rectPair2) == rectPair4)

testReport testName result = putStrLn $ (if result then "Yes " else "No ") ++ testName