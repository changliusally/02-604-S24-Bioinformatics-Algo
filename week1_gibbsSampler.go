// 1.2.9 GibbsSampler (GO version)
package main
import (
    "fmt"
    "bufio"
    "os"
    "strings"
    "strconv"
    "math/rand" //This should be helpful
  )

//Insert your GibbsSampler function here, along with any subroutines you need
func GibbsSampler(Dna []string, k, t, N int) []string {
    bestscore := k*t
    bestmotif := make([]string, t)
    // generate random initial motifs 20 times, keep the motifs with the smallest score
    for i := 0; i < 20; i++ {
        motifs := GibbsSamplerOnce(Dna, k, t, N)
        if score(motifs) < bestscore {
            bestscore = score(motifs)
            bestmotif = motifs
        }
    }
    return bestmotif
}

// run a single GibbsSampler with only one random initial motifs, and do N times of random selection of updated motifs
func GibbsSamplerOnce(Dna []string, k, t, N int) []string {
    initial := RandomizedInitial(Dna, k, t)
    best := initial

    for j:=0; j<N; j++ {
        // choose one of the t motifs randomly and ignored when calculating profile matrix
        r := rand.Intn(t)
        motif_profile := make([]string, t-1)
        // create a list of dna except for the chosen motif r
        if r != 0 && r != t-1 {
            motif_profile = append(best[:r],best[r+1:]...)
        } else if r == 0 {
            motif_profile = best[1:]
        } else {
            motif_profile = best[:t-1]
        }
        // calculate the profile matrix of the list of dna excluding motif r
        profile := Profile(motif_profile)
 		// based on the preceding profile matrix, do a weighted random selection of the k-mers of motif r
        motif := RandomMotif(profile, Dna[r])
		// put the motif back into the right position of the motif list (in the orginal order)
        Motif := make([]string, t)
        // if the removed DNA seq is not at the start and end, then we can insert it by append the motif ahead and behind
        if r != 0 && r != t-1 {
            Motif = append(best[:r],motif)
            Motif = append(Motif, best[r+1:]...)
        } else if r == 0 { // if it is the first motif
            Motif = make([]string, 0)
            Motif = append(Motif, motif)
            Motif = append(Motif, motif_profile...)
        } else { // if it is the last motif
            for k := range motif_profile {
                Motif[k] = motif_profile[k]
            }
            Motif[r] = motif
        }
		// check if the updated motif is improving
        if score(Motif) < score(best){
            best = Motif
        }
    }
    return best
}

// calculate the probability of each k-mer of motif r
func calProbability(profile [][]float64, text string) []float64 {
    probability := make([]float64,0)
    // traverse every possible start point
    for i := 0; i < len(text)-len(profile[0])+1; i++ {
        prob := 1.0
        for j := 0; j < len(profile[0]); j++ {
            if text[i+j] == 'A' {
    prob = prob * profile[0][j]
   } else if text[i+j] == 'C' {
    prob = prob * profile[1][j]
   } else if text[i+j] == 'G' {
    prob = prob * profile[2][j] 
   } else {
    prob = prob * profile[3][j]
   }
        }
        probability = append(probability, prob)
    }
    return probability
}


// based on the probability matrix, do a weighted random selection of k-mer of motif r
func RandomMotif(profile [][]float64, text string) string {
    // generate the weight matrix
    prob := calProbability(profile, text)
 sum := 0.0
    cumulativeProbabilities := make([]float64, len(text)-len(profile[0])+1)
 for i, weight := range prob {
  sum += weight
  cumulativeProbabilities[i] = sum
 }


 // Generate a random number in the range [0, sum)
 randomNumber := rand.Float64() * sum

 // Find the index of the first cumulative probability greater than the random number
 for i, cumulativeProbability := range cumulativeProbabilities {
  if randomNumber < cumulativeProbability {
            return text[i:i+len(profile[0])]
  }
 }
    return text[0:1]
}


// generate a random k-mer for each DNA
func RandomizedInitial(Dna []string, k int, t int) []string {
    // motifs is a t*k string matrix 
    initial := make([]string, 0)
    for i:= 0; i < t; i++{
        index := rand.Intn(len(Dna[i])-k+1)
        motif := Dna[i][index:index+k]
  initial = append(initial, motif)
    }
 return initial
}


// calculate the score of a give motif
func score(motifs []string) int {

	score := 0
	// for each column
	for i := range motifs[0] {

		Acount := 0
		Ccount := 0
		Gcount := 0
		Tcount := 0
		// for each row
		for j := range motifs {
		if motifs[j][i] == 'A' {
			Acount += 1
		} else if motifs[j][i] == 'C' {
			Ccount += 1
		} else if motifs[j][i] == 'G' {
			Gcount += 1
		} else {
			Tcount += 1
		}
		}
		max := findMax(Acount, Ccount, Gcount, Tcount)
		diff := Acount + Ccount + Gcount + Tcount - max
		score += diff
 	}
 	return score
}


func findMax(a, b, c, d int) int {
 max := a

 if b > max {
  max = b
 }
 if c > max {
  max = c
 }
 if d > max {
  max = d
 }

 return max
}


// calculate the profile matrix for the motif set
func Profile(motifs []string) [][]float64 {
    // build a 4*k matrix to store probability
    profile := make([][]float64, 4)
 for k := range profile {
  profile[k] = make([]float64, len(motifs[0])) 
 }
    // for each column
 for i := range motifs[0]{
  // considering pseudocount
  Acount := 1
  Ccount := 1
  Gcount := 1
  Tcount := 1
  // for each row
  for j := range motifs {
   if motifs[j][i] == 'A' {
    Acount += 1
   } else if motifs[j][i] == 'C' {
    Ccount += 1
   } else if motifs[j][i] == 'G' {
    Gcount += 1
   } else {
    Tcount += 1
   }
  }
        // total := len(motifs)
        total := float64(Acount + Ccount + Gcount + Tcount)
  
  profile[0][i] = float64(Acount)/total
  profile[1][i] = float64(Ccount)/total
  profile[2][i] = float64(Gcount)/total
  profile[3][i] = float64(Tcount)/total
  
 }
 return profile
}