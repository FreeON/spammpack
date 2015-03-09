import MatrixTree
import System.CPUTime
import System.Random (getStdGen, randomRs, StdGen)

main = do let sizes = fmap (2^) [1..20]
          putStrLn "norm"
          mapM_ (doTiming setNorm) sizes

doTiming :: (MatrixTree -> Norm) -> Int -> IO ()
doTiming op size = do gen <- getStdGen ; let tree = makeRandomTree size gen
                      t1 <- getCPUTime
                      norm <- return $ op tree
                      t2 <- getCPUTime
                      let time = fromIntegral (t2 - t1) * 1e-12 :: Double
                      putStrLn $ show norm ++ " " ++ show size ++ " " ++ show time

makeRandomTree :: Int -> StdGen -> MatrixTree
makeRandomTree size gen = indexedListToTree (size, size, ijxs)
                          where ijxs = zipWith (\(i, j) x -> (i, j, x)) indices randomNums
                                indices = [(i, j) | j <- [1..size], i <- [1..size]]
                                randomNums = take (size^2) $ randomRs (0.0, 1.0) gen
