import Data.List (intersperse)
import Data.Time.Clock (getCurrentTime, diffUTCTime)
import MatrixTree (matrixListToTree, MatrixTree, MTree, Norm, norm)
import SpAMM (treeAdd, treeMult)
import System.Directory (removeFile)
import System.IO (Handle, hClose, hPutStr, openTempFile)
import System.Random (getStdGen, newStdGen, randomRs, StdGen)

main = do (tempName, tempHandle) <- openTempFile "." "temp"
          let sizes = fmap (2^) [5..10]
          putStrLn "treeAdd"
          mapM_ (doTiming tempHandle treeAdd) sizes
          putStrLn "treeMult"
          mapM_ (doTiming tempHandle treeMult) sizes
          hClose tempHandle; removeFile tempName

doTiming :: Handle -> (MatrixTree -> MatrixTree -> MatrixTree) -> Int -> IO ()
doTiming handle op size =
         do gen <- getStdGen ; let fstTree = makeRandomTree size gen
            gen <- newStdGen ; let sndTree = makeRandomTree size gen
            let (h1, w1, mTree1) = fstTree ; (h2, w2, mTree2) = sndTree
            hPutStr handle $ show (norm mTree1) ; hPutStr handle $ show (norm mTree2)
            t1 <- getCurrentTime
            let tree = op fstTree sndTree
            let (h, w, mTree) = tree
            hPutStr handle (show $ norm mTree)
            t2 <- getCurrentTime
            printList [show size, show $ diffUTCTime t2 t1]

makeRandomTree :: Int -> StdGen -> MatrixTree
makeRandomTree size gen = matrixListToTree (size, size, ijxs)
               where ijxs = zipWith (\(i, j) x -> (i, j, x)) indices randomNums
                     indices = [(i, j) | j <- [1..size], i <- [1..size]]
                     randomNums = take (size^2) $ randomRs (0.0, 1.0) gen

printList :: [String] -> IO ()
printList = putStrLn . concat . intersperse " "