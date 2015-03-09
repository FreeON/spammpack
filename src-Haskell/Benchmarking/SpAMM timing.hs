import MatrixTree ; import SpAMM ; import System.CPUTime ; import System.Directory (removeFile)
import System.IO (hClose, openTempFile)

main = do let sizes = fmap (2^) [11..20]
          putStrLn "treeAdd"
          mapM_ (doTiming treeAdd) sizes
          putStrLn "treeMult"
          mapM_ (doTiming treeMult) sizes

doTiming :: (MatrixTree -> MatrixTree -> MatrixTree) -> Int -> IO ()
doTiming op size = do let fstTree = makeRandomTree size ; sndTree = makeRandomTree size
                      t1 <- getCPUTime
                      tree <- return $ op fstTree sndTree
                      t2 <- getCPUTime
                      (tempName, tempHandle) <- openTempFile "." "temp"
                      hClose tempHandle
                      removeFile tempName
                      writeTreeToMatrixMarket tree "coordinate" tempName
                      removeFile tempName
                      let time = fromIntegral (t2 - t1) * 1e-12 :: Double
                      putStrLn $ show size ++ " " ++ show time

makeRandomTree :: Int -> MatrixTree
makeRandomTree size = indexedListToTree (size, size, ijxs)
                      where ijxs = [(i, j, x) | j <- [1..size], i <- [1..size], x <- randomNums]
                            randomNums = take (size^2) $ randomRs (0.0, 1.0) (mkStdGen 245)
