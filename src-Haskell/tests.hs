import SpAMM

main = do testRowsSparse ; testAdd ; testMult

testRowsSparse = do
                 testRowListTree <- readTreeFromRowList "testrowlist.txt" 
                 testSparseListTree <- readTreeFromIndexedList "testindexedlist.txt" 
                 print $ testRowListTree == testSparseListTree

testAdd = do
          testAddTree <- readTreeFromRowList "testaddlist.txt"
          testAddResultTree <- readTreeFromRowList "testaddresultlist.txt"
          print $ testAddTree `treeAdd` testAddTree == testAddResultTree

testMult = do
           testMultTree1 <- readTreeFromRowList "testmultlist1.txt"
           testMultTree2 <- readTreeFromRowList "testmultlist2.txt"
           testMultResultTree <- readTreeFromRowList "testmultresultlist.txt"
           print $ testMultTree1 `treeMult` testMultTree2 == testMultResultTree