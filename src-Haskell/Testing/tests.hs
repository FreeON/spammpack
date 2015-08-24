import MatrixTree ; import SpAMM ; import System.Directory (removeFile)

main = do testReadWrite ; testAdd ; testMultiply

testReadWrite = do tree1 <- mmReadTree "Matrices/rect4arr.txt"
                   mmWriteTree tree1 "coordinate" "Matrices/tempfile.txt"
                   tree2 <- mmReadTree "Matrices/tempfile.txt"
                   testReport "readWrite" (tree1 == tree2)
                   removeFile "Matrices/tempfile.txt"

testAdd = do addZeros ; addZeroToValue ; addZeroToRect;
             addValueToValue ; addRectToRect

addZeros = do zeroTree <- mmReadTree "Matrices/zerorectcoor.txt"
              testReport "addZeros" (zeroTree `treeAdd` zeroTree == zeroTree)

addZeroToValue = do valueTree <- mmReadTree "Matrices/value1coor.txt"
                    zeroTree <- mmReadTree "Matrices/zerovaluecoor.txt"
                    testReport "addValuetoZero" (valueTree `treeAdd` zeroTree
                                                 == valueTree)
                    testReport "addZerotoValue" (zeroTree `treeAdd` valueTree
                                                 == valueTree)

addZeroToRect = do rectTree <- mmReadTree "Matrices/rect1coor.txt"
                   zeroTree <- mmReadTree "Matrices/zerorectcoor.txt"
                   testReport "addRectToZero" (rectTree `treeAdd` zeroTree
                                               == rectTree)
                   testReport "addZeroToRect" (zeroTree `treeAdd` rectTree
                                               == rectTree)

addValueToValue = do valueTree1 <- mmReadTree "Matrices/value1coor.txt"
                     valueTree2 <- mmReadTree "Matrices/value2coor.txt"
                     valueTree3 <- mmReadTree "Matrices/value3coor.txt"
                     testReport "addValueToValue" (valueTree1 `treeAdd`
                                                   valueTree2 ==
                                                   valueTree3)

addRectToRect = do rectTree1 <- mmReadTree "Matrices/rect1coor.txt"
                   rectTree2 <- mmReadTree "Matrices/rect2coor.txt"
                   rectTree3 <- mmReadTree "Matrices/rect3coor.txt"
                   testReport "addRectToRect" (rectTree1 `treeAdd` rectTree2
                                               == rectTree3)

testMultiply = do multZeros ; multZeroByValue ; multZeroByRect ;
                  multValueByValue ; multValueByRect ; multRectByRect ;
                  multRectByRectArray

multZeros = do zeroTree <- mmReadTree "Matrices/zerorectcoor.txt"
               zeroSquare <- mmReadTree "Matrices/zerosquarecoor.txt"
               testReport "multZeros" (zeroTree `treeMult`
                                       (treeTranspose zeroTree)
                                       == zeroSquare)

multZeroByValue = do zeroRow <- mmReadTree "Matrices/zerorowcoor.txt"
                     valueTree <- mmReadTree "Matrices/value1coor.txt"
                     testReport "multValueByZero" (valueTree `treeMult`
                                                   zeroRow == zeroRow)
                     testReport "multZeroByValue" ((treeTranspose zeroRow)
                                                   `treeMult`
                                                   (treeTranspose valueTree)
                                                   == treeTranspose zeroRow)

multZeroByRect = do zeroTree <- mmReadTree "Matrices/zerorectcoor.txt"
                    rectTree <- mmReadTree "Matrices/rect1coor.txt"
                    zeroSquare <- mmReadTree "Matrices/zerosquarecoor.txt"
                    testReport "multZeroByRect" (zeroTree `treeMult`
                                                 (treeTranspose rectTree)
                                                 == zeroSquare)
                    testReport "multZeroByRect" (rectTree `treeMult`
                                                 (treeTranspose zeroTree)
                                                  == zeroSquare)

multValueByValue =
    do valueTree1 <- mmReadTree "Matrices/value2coor.txt"
       valueTree2 <- mmReadTree "Matrices/value3coor.txt"
       valueTree3 <- mmReadTree "Matrices/value4coor.txt"
       testReport "multValueByValue" (valueTree1 `treeMult` valueTree2
                                      == valueTree3)

multValueByRect =
    do valueTree <- mmReadTree "Matrices/value2coor.txt"
       rowTree <- mmReadTree "Matrices/rowcoor.txt"
       valuetimesrowTree <- mmReadTree "Matrices/valuetimesrowcoor.txt"
       testReport "multValueByRect" (valueTree `treeMult` rowTree
                                     == valuetimesrowTree)
       testReport "multRectByValue" ((treeTranspose rowTree) `treeMult`
                                     (treeTranspose valueTree) ==
                                     treeTranspose valuetimesrowTree)

multRectByRect =
    do rectTree1 <- mmReadTree "Matrices/rect1coor.txt"
       rectTree2 <- mmReadTree "Matrices/rect2coor.txt"
       rectTree4 <- mmReadTree "Matrices/rect4coor.txt"
       testReport "multRectByRect" (rectTree1 `treeMult`
                                    (treeTranspose rectTree2) == rectTree4)

multRectByRectArray =
    do rectTree1 <- mmReadTree "Matrices/rect1arr.txt"
       rectTree2 <- mmReadTree "Matrices/rect2arr.txt"
       rectTree4 <- mmReadTree "Matrices/rect4arr.txt"
       testReport "multRectByRectArray" (rectTree1 `treeMult`
                                         (treeTranspose rectTree2)
                                         == rectTree4)

testReport testName result =
           putStrLn $ (if result then "Yes " else "No ") ++ testName