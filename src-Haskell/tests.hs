import SpAMM

main = do testAdd ; testMultiply

testAdd = do addZeros ; addZeroToValue ; addZeroToRect;
             addValueToValue ; addRectToRect

addZeros = do zeroTree <- readTreeFromRowList "Test matrices/zerorect.txt"
              testReport "addZeros" (zeroTree `treeAdd` zeroTree == zeroTree)

addZeroToValue = do valueTree <- readTreeFromRowList "Test matrices/value.txt"
                    zeroTree <- readTreeFromRowList "Test matrices/zerovalue.txt"
                    testReport "addValuetoZero" (valueTree `treeAdd` zeroTree == valueTree)
                    testReport "addZerotoValue" (zeroTree `treeAdd` valueTree == valueTree)

addZeroToRect = do rectTree <- readTreeFromRowList "Test matrices/rect.txt"
                   zeroTree <- readTreeFromRowList "Test matrices/zerorect.txt"
                   testReport "addRectToZero" (rectTree `treeAdd` zeroTree == rectTree)
                   testReport "addZeroToRect" (zeroTree `treeAdd` rectTree == rectTree)

addValueToValue = do valueTree1 <- readTreeFromRowList "Test matrices/value.txt"
                     valueTree2 <- readTreeFromRowList "Test matrices/value2.txt"
                     valueTree3 <- readTreeFromRowList "Test matrices/value3.txt"
                     testReport "addValueToValue" (valueTree1 `treeAdd` valueTree2 == valueTree3)

addRectToRect = do rectTree1 <- readTreeFromRowList "Test matrices/rect.txt"
                   rectTree2 <- readTreeFromRowList "Test matrices/rect2.txt"
                   rectTree3 <- readTreeFromRowList "Test matrices/rect3.txt"
                   testReport "addRectToRect" (rectTree1 `treeAdd` rectTree2 == rectTree3)

testMultiply = do multZeros ; multZeroByValue ; multZeroByRect ;
                  multValueByValue ; multValueByRect ; multRectByRect

multZeros = do zeroTree <- readTreeFromRowList "Test matrices/zerorect.txt"
               zeroSquare <- readTreeFromRowList "Test matrices/zerosquare.txt"
               testReport "multZeros" (zeroTree `treeMult` (treeTranspose zeroTree) == zeroSquare)

multZeroByValue = do zeroRow <- readTreeFromRowList "Test matrices/zerorow.txt"
                     valueTree <- readTreeFromRowList "Test matrices/value.txt"
                     testReport "multValueByZero" (valueTree `treeMult` zeroRow == zeroRow)
                     testReport "multZeroByValue" ((treeTranspose zeroRow) `treeMult`
                                                  (treeTranspose valueTree) == treeTranspose zeroRow)

multZeroByRect = do zeroTree <- readTreeFromRowList "Test matrices/zerorect.txt"
                    rectTree <- readTreeFromRowList "Test matrices/rect.txt"
                    zeroSquare <- readTreeFromRowList "Test matrices/zerosquare.txt"
                    testReport "multZeroByRect" (zeroTree `treeMult`
                                                (treeTranspose rectTree) == zeroSquare)
                    testReport "multZeroByRect" (rectTree `treeMult` (treeTranspose zeroTree)
                                                 == zeroSquare)

multValueByValue = do valueTree1 <- readTreeFromRowList "Test matrices/value2.txt"
                      valueTree2 <- readTreeFromRowList "Test matrices/value3.txt"
                      valueTree3 <- readTreeFromRowList "Test matrices/value6.txt"
                      testReport "multValueByValue" (valueTree1 `treeMult` valueTree2 == valueTree3)

multValueByRect = do valueTree <- readTreeFromRowList "Test matrices/value2.txt"
                     rowTree <- readTreeFromRowList "Test matrices/row.txt"
                     valuetimesrowTree <- readTreeFromRowList "Test matrices/valuetimesrow.txt"
                     testReport "multValueByRect" (valueTree `treeMult` rowTree == valuetimesrowTree)
                     testReport "multRectByValue" ((treeTranspose rowTree) `treeMult`
                                                   (treeTranspose valueTree) ==
                                                    treeTranspose valuetimesrowTree)

multRectByRect = do rectTree1 <- readTreeFromRowList "Test matrices/rect.txt"
                    rectTree2 <- readTreeFromRowList "Test matrices/rect2.txt"
                    rectTree4 <- readTreeFromRowList "Test matrices/rect4.txt"
                    testReport "multRectByRect" (rectTree1 `treeMult`
                                                (treeTranspose rectTree2) == rectTree4)

testReport testName result = putStrLn $ (if result then "Yes " else "No ") ++ testName