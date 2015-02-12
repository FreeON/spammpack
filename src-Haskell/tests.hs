import MatrixMarket ; import SpAMM

main = do testAdd ; testMultiply

testAdd = do addZeros ; addZeroToValue ; addZeroToRect;
             addValueToValue ; addRectToRect

addZeros = do zeroTree <- readTreeFromMatrixMarket "Test matrices/zerorectcoor.txt"
              testReport "addZeros" (zeroTree `treeAdd` zeroTree == zeroTree)

addZeroToValue = do valueTree <- readTreeFromMatrixMarket "Test matrices/value1coor.txt"
                    zeroTree <- readTreeFromMatrixMarket "Test matrices/zerovaluecoor.txt"
                    testReport "addValuetoZero" (valueTree `treeAdd` zeroTree == valueTree)
                    testReport "addZerotoValue" (zeroTree `treeAdd` valueTree == valueTree)

addZeroToRect = do rectTree <- readTreeFromMatrixMarket "Test matrices/rect1coor.txt"
                   zeroTree <- readTreeFromMatrixMarket "Test matrices/zerorectcoor.txt"
                   testReport "addRectToZero" (rectTree `treeAdd` zeroTree == rectTree)
                   testReport "addZeroToRect" (zeroTree `treeAdd` rectTree == rectTree)

addValueToValue = do valueTree1 <- readTreeFromMatrixMarket "Test matrices/value1coor.txt"
                     valueTree2 <- readTreeFromMatrixMarket "Test matrices/value2coor.txt"
                     valueTree3 <- readTreeFromMatrixMarket "Test matrices/value3coor.txt"
                     testReport "addValueToValue" (valueTree1 `treeAdd` valueTree2 == valueTree3)

addRectToRect = do rectTree1 <- readTreeFromMatrixMarket "Test matrices/rect1coor.txt"
                   rectTree2 <- readTreeFromMatrixMarket "Test matrices/rect2coor.txt"
                   rectTree3 <- readTreeFromMatrixMarket "Test matrices/rect3coor.txt"
                   testReport "addRectToRect" (rectTree1 `treeAdd` rectTree2 == rectTree3)

testMultiply = do multZeros ; multZeroByValue ; multZeroByRect ; multValueByValue ;
                  multValueByRect ; multRectByRect ; multRectByRectArray

multZeros = do zeroTree <- readTreeFromMatrixMarket "Test matrices/zerorectcoor.txt"
               zeroSquare <- readTreeFromMatrixMarket "Test matrices/zerosquarecoor.txt"
               testReport "multZeros" (zeroTree `treeMult` (treeTranspose zeroTree) == zeroSquare)

multZeroByValue = do zeroRow <- readTreeFromMatrixMarket "Test matrices/zerorowcoor.txt"
                     valueTree <- readTreeFromMatrixMarket "Test matrices/value1coor.txt"
                     testReport "multValueByZero" (valueTree `treeMult` zeroRow == zeroRow)
                     testReport "multZeroByValue" ((treeTranspose zeroRow) `treeMult`
                                                  (treeTranspose valueTree) == treeTranspose zeroRow)

multZeroByRect = do zeroTree <- readTreeFromMatrixMarket "Test matrices/zerorectcoor.txt"
                    rectTree <- readTreeFromMatrixMarket "Test matrices/rect1coor.txt"
                    zeroSquare <- readTreeFromMatrixMarket "Test matrices/zerosquarecoor.txt"
                    testReport "multZeroByRect" (zeroTree `treeMult`
                                                (treeTranspose rectTree) == zeroSquare)
                    testReport "multZeroByRect" (rectTree `treeMult` (treeTranspose zeroTree)
                                                 == zeroSquare)

multValueByValue = do valueTree1 <- readTreeFromMatrixMarket "Test matrices/value2coor.txt"
                      valueTree2 <- readTreeFromMatrixMarket "Test matrices/value3coor.txt"
                      valueTree3 <- readTreeFromMatrixMarket "Test matrices/value4coor.txt"
                      testReport "multValueByValue" (valueTree1 `treeMult` valueTree2 == valueTree3)

multValueByRect = do valueTree <- readTreeFromMatrixMarket "Test matrices/value2coor.txt"
                     rowTree <- readTreeFromMatrixMarket "Test matrices/rowcoor.txt"
                     valuetimesrowTree <- readTreeFromMatrixMarket "Test matrices/valuetimesrowcoor.txt"
                     testReport "multValueByRect" (valueTree `treeMult` rowTree == valuetimesrowTree)
                     testReport "multRectByValue" ((treeTranspose rowTree) `treeMult`
                                                   (treeTranspose valueTree) ==
                                                    treeTranspose valuetimesrowTree)

multRectByRect = do rectTree1 <- readTreeFromMatrixMarket "Test matrices/rect1coor.txt"
                    rectTree2 <- readTreeFromMatrixMarket "Test matrices/rect2coor.txt"
                    rectTree4 <- readTreeFromMatrixMarket "Test matrices/rect4coor.txt"
                    testReport "multRectByRect" (rectTree1 `treeMult`
                                                (treeTranspose rectTree2) == rectTree4)

multRectByRectArray = do rectTree1 <- readTreeFromMatrixMarket "Test matrices/rect1arr.txt"
                         rectTree2 <- readTreeFromMatrixMarket "Test matrices/rect2arr.txt"
                         rectTree4 <- readTreeFromMatrixMarket "Test matrices/rect4arr.txt"
                         testReport "multRectByRectArray" (rectTree1 `treeMult`
                                                          (treeTranspose rectTree2) == rectTree4)

testReport testName result = putStrLn $ (if result then "Yes " else "No ") ++ testName