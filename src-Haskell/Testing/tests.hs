import MatrixTree ; import SpAMM ; import System.Directory (removeFile)

main = do testReadWrite ; testAdd ; testMultiply

testReadWrite = do tree1 <- readTreeFromMatrixMarket "Matrices/rect4arr.txt"
                   writeTreeToMatrixMarket tree1 "coordinate" "Matrices/tempfile.txt"
                   tree2 <- readTreeFromMatrixMarket "Matrices/tempfile.txt"
                   testReport "readWrite" (tree1 == tree2)
                   removeFile "Matrices/tempfile.txt"

testAdd = do addZeros ; addZeroToValue ; addZeroToRect;
             addValueToValue ; addRectToRect

addZeros = do zeroTree <- readTreeFromMatrixMarket "Matrices/zerorectcoor.txt"
              testReport "addZeros" (zeroTree `treeAdd` zeroTree == zeroTree)

addZeroToValue = do valueTree <- readTreeFromMatrixMarket "Matrices/value1coor.txt"
                    zeroTree <- readTreeFromMatrixMarket "Matrices/zerovaluecoor.txt"
                    testReport "addValuetoZero" (valueTree `treeAdd` zeroTree == valueTree)
                    testReport "addZerotoValue" (zeroTree `treeAdd` valueTree == valueTree)

addZeroToRect = do rectTree <- readTreeFromMatrixMarket "Matrices/rect1coor.txt"
                   zeroTree <- readTreeFromMatrixMarket "Matrices/zerorectcoor.txt"
                   testReport "addRectToZero" (rectTree `treeAdd` zeroTree == rectTree)
                   testReport "addZeroToRect" (zeroTree `treeAdd` rectTree == rectTree)

addValueToValue = do valueTree1 <- readTreeFromMatrixMarket "Matrices/value1coor.txt"
                     valueTree2 <- readTreeFromMatrixMarket "Matrices/value2coor.txt"
                     valueTree3 <- readTreeFromMatrixMarket "Matrices/value3coor.txt"
                     testReport "addValueToValue" (valueTree1 `treeAdd` valueTree2 == valueTree3)

addRectToRect = do rectTree1 <- readTreeFromMatrixMarket "Matrices/rect1coor.txt"
                   rectTree2 <- readTreeFromMatrixMarket "Matrices/rect2coor.txt"
                   rectTree3 <- readTreeFromMatrixMarket "Matrices/rect3coor.txt"
                   testReport "addRectToRect" (rectTree1 `treeAdd` rectTree2 == rectTree3)

testMultiply = do multZeros ; multZeroByValue ; multZeroByRect ; multValueByValue ;
                  multValueByRect ; multRectByRect ; multRectByRectArray

multZeros = do zeroTree <- readTreeFromMatrixMarket "Matrices/zerorectcoor.txt"
               zeroSquare <- readTreeFromMatrixMarket "Matrices/zerosquarecoor.txt"
               testReport "multZeros" (zeroTree `treeMult` (treeTranspose zeroTree) == zeroSquare)

multZeroByValue = do zeroRow <- readTreeFromMatrixMarket "Matrices/zerorowcoor.txt"
                     valueTree <- readTreeFromMatrixMarket "Matrices/value1coor.txt"
                     testReport "multValueByZero" (valueTree `treeMult` zeroRow == zeroRow)
                     testReport "multZeroByValue" ((treeTranspose zeroRow) `treeMult`
                                                  (treeTranspose valueTree) == treeTranspose zeroRow)

multZeroByRect = do zeroTree <- readTreeFromMatrixMarket "Matrices/zerorectcoor.txt"
                    rectTree <- readTreeFromMatrixMarket "Matrices/rect1coor.txt"
                    zeroSquare <- readTreeFromMatrixMarket "Matrices/zerosquarecoor.txt"
                    testReport "multZeroByRect" (zeroTree `treeMult`
                                                (treeTranspose rectTree) == zeroSquare)
                    testReport "multZeroByRect" (rectTree `treeMult` (treeTranspose zeroTree)
                                                 == zeroSquare)

multValueByValue = do valueTree1 <- readTreeFromMatrixMarket "Matrices/value2coor.txt"
                      valueTree2 <- readTreeFromMatrixMarket "Matrices/value3coor.txt"
                      valueTree3 <- readTreeFromMatrixMarket "Matrices/value4coor.txt"
                      testReport "multValueByValue" (valueTree1 `treeMult` valueTree2 == valueTree3)

multValueByRect = do valueTree <- readTreeFromMatrixMarket "Matrices/value2coor.txt"
                     rowTree <- readTreeFromMatrixMarket "Matrices/rowcoor.txt"
                     valuetimesrowTree <- readTreeFromMatrixMarket "Matrices/valuetimesrowcoor.txt"
                     testReport "multValueByRect" (valueTree `treeMult` rowTree == valuetimesrowTree)
                     testReport "multRectByValue" ((treeTranspose rowTree) `treeMult`
                                                   (treeTranspose valueTree) ==
                                                    treeTranspose valuetimesrowTree)

multRectByRect = do rectTree1 <- readTreeFromMatrixMarket "Matrices/rect1coor.txt"
                    rectTree2 <- readTreeFromMatrixMarket "Matrices/rect2coor.txt"
                    rectTree4 <- readTreeFromMatrixMarket "Matrices/rect4coor.txt"
                    testReport "multRectByRect" (rectTree1 `treeMult`
                                                (treeTranspose rectTree2) == rectTree4)

multRectByRectArray = do rectTree1 <- readTreeFromMatrixMarket "Matrices/rect1arr.txt"
                         rectTree2 <- readTreeFromMatrixMarket "Matrices/rect2arr.txt"
                         rectTree4 <- readTreeFromMatrixMarket "Matrices/rect4arr.txt"
                         testReport "multRectByRectArray" (rectTree1 `treeMult`
                                                          (treeTranspose rectTree2) == rectTree4)

testReport testName result = putStrLn $ (if result then "Yes " else "No ") ++ testName