import MatrixTree ; import System.CPUTime ; import System.Random (mkStdGen, randomRs)

main = do let sizes = fmap (2^) [11..20]
          putStrLn "norm"
          mapM_ (doTiming setNorm) sizes

doTiming :: (MatrixTree -> Norm) -> Int -> IO ()
doTiming op size = do let tree = makeRandomTree size
                      t1 <- getCPUTime
                      norm <- return $ op tree
                      t2 <- getCPUTime
                      let time = fromIntegral (t2 - t1) * 1e-12 :: Double
                      putStrLn $ show norm ++ " " ++ show size ++ " " ++ show time

makeRandomTree :: Int -> MatrixTree
makeRandomTree size = indexedListToTree (size, size, ijxs)
                      where ijxs = [(i, j, x) | j <- [1..size], i <- [1..size], x <- randomNums]
                            randomNums = take (size^2) $ randomRs (0.0, 1.0) (mkStdGen 245)
