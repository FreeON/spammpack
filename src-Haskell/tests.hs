import SpAMM

main = do testRowsSparse ; testAdd ; testMult

testRowsSparse = do
                 testRowListTree <- readListToTree "testrowlist.txt" 
                 testSparseListTree <- readIndexedListToTree "testindexedlist.txt" 
                 print $ testRowListTree == testSparseListTree

testAdd = do
          testAddTree <- readListToTree "testaddlist.txt"
          testAddResultTree <- readListToTree "testaddresultlist.txt"
          print $ testAddTree `treeAdd` testAddTree == testAddResultTree

testMult = do
           testMultTree1 <- readListToTree "testmultlist1.txt"
           testMultTree2 <- readListToTree "testmultlist2.txt"
           testMultResultTree <- readListToTree "testmultresultlist.txt"
           print $ testMultTree1 `treeMult` testMultTree2 == testMultResultTree

readListToTree filePath = do list <- readMatrixList filePath
                             return $ listToTree list

readIndexedListToTree filePath = do indexedList <- readIndexedList filePath
                                    return $ indexedListToTree indexedList