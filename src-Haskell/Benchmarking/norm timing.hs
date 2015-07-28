import Data.List (intersperse)
import Data.Time.Clock (getCurrentTime, diffUTCTime)
import MatrixTree (matrixListToTree, MatrixTree, MTree, Norm, norm, setNorm)
import System.Directory (removeFile)
import System.IO (Handle, hClose, hPutStr, openTempFile)
import System.Random (getStdGen, randomRs, StdGen)

main = do (tempName, tempHandle) <- openTempFile "." "temp"
          let sizes = fmap (2^) [5..11]
          putStrLn "norm"
          mapM_ (doTiming tempHandle (norm . setNorm)) sizes
          hClose tempHandle; removeFile tempName

doTiming :: Handle -> (MTree -> Norm) -> Int -> IO ()
doTiming handle op size =
         do gen <- getStdGen ; let (h, w, tree) = makeRandomTree size gen
            hPutStr handle $ show (norm tree)
            t1 <- getCurrentTime
            let num = op tree
            hPutStr handle $ show num
            t2 <- getCurrentTime
            printList [show size, show $ diffUTCTime t2 t1]

makeRandomTree :: Int -> StdGen -> MatrixTree
makeRandomTree size gen = matrixListToTree (size, size, ijxs)
               where ijxs = zipWith (\(i, j) x -> (i, j, x)) indices randomNums
                     indices = [(i, j) | j <- [1..size], i <- [1..size]]
                     randomNums = take (size^2) $ randomRs (0.0, 1.0) gen

printList :: [String] -> IO ()
printList = putStrLn . concat . intersperse " "