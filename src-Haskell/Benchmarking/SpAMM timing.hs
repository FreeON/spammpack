import MatrixTree ; import SpAMM
import System.CPUTime
import System.Random (getStdGen, newStdGen, randomRs, StdGen)

main = do let sizes = fmap (2^) [1..20]
          putStrLn "treeAdd"
          mapM_ (doTiming treeAdd) sizes
          putStrLn "treeMult"
          mapM_ (doTiming treeMult) sizes

doTiming :: (MatrixTree -> MatrixTree -> MatrixTree) -> Int -> IO ()
doTiming op size = do gen <- getStdGen ; let fstTree = makeRandomTree size gen
                      gen <- newStdGen ; let sndTree = makeRandomTree size gen
                      t1 <- getCPUTime
                      tree <- return $ op fstTree sndTree
                      t2 <- getCPUTime
                      let time = fromIntegral (t2 - t1) * 1e-12 :: Double
                      putStrLn $ show (norm tree) ++ " " ++ show size ++ " " ++ show time

makeRandomTree :: Int -> StdGen -> MatrixTree
makeRandomTree size gen = indexedListToTree (size, size, ijxs)
                          where ijxs = zipWith (\(i, j) x -> (i, j, x)) indices randomNums
                                indices = [(i, j) | j <- [1..size], i <- [1..size]]
                                randomNums = take (size^2) $ randomRs (0.0, 1.0) gen
