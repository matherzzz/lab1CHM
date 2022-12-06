package main

import (
	"fmt"
	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/opts"
	"github.com/go-echarts/go-echarts/v2/types"
	"math"
	"math/rand"
	"os"
)

const (
	a = -1.0
	b = 1.0
	n = 10
	m = n + 2
)

func generateLineItems(rhoVec []float64) []opts.LineData {
	items := make([]opts.LineData, 0)
	for i := 0; i < len(rhoVec); i++ {
		items = append(items, opts.LineData{Value: rhoVec[i]})
	}
	return items
}

func createLineChart(rhoVec []float64) {
	// create a new line instance
	line := charts.NewLine()
	// set some global options like Title/Legend/ToolTip or anything else
	line.SetGlobalOptions(
		charts.WithInitializationOpts(opts.Initialization{
			Theme: types.ThemeInfographic}),
		charts.WithTitleOpts(opts.Title{
			Title:    "Lab1 in Go",
			Subtitle: "pho",
		}))
	// Put data into instance
	t := make([]int, 0)
	for i := 0; i < len(rhoVec); i++ {
		t = append(t, i)
	}
	line.SetXAxis(t).
		AddSeries("Category A", generateLineItems(rhoVec)).
		SetSeriesOptions(charts.WithLineChartOpts(opts.LineChart{Smooth: true}))
	f, _ := os.Create("rho.html")
	_ = line.Render(f)
}

func normMtr(A [m][n]float64) float64 {
	var norm float64
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			norm += A[i][j] * A[i][j]
		}
	}
	return math.Sqrt(norm)
}

func normVecM(f [m]float64) float64 {
	var norm float64
	for i := 0; i < m; i++ {
		norm += f[i] * f[i]
	}
	return math.Sqrt(norm)
}

func normVecN(f [n]float64) float64 {
	var norm float64
	for i := 0; i < n; i++ {
		norm += f[i] * f[i]
	}
	return math.Sqrt(norm)
}

func mulMtrVecMN(A [m][n]float64, u [n]float64) [m]float64 {
	var f [m]float64
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			f[i] += A[i][j] * u[j]
		}
	}
	return f
}

func mulVecVecNN(f [n]float64, u [n]float64) float64 {
	var answer float64
	for i := 0; i < n; i++ {
		answer += f[i] * u[i]
	}
	return answer
}

func mulMtrVecNM(A [n][m]float64, f [m]float64) [n]float64 {
	var answer [n]float64
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			answer[i] += A[i][j] * f[j]
		}
	}
	return answer
}

func mulMtrVecNN(A [n][n]float64, f [n]float64) [n]float64 {
	var answer [n]float64
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			answer[i] += A[i][j] * f[j]
		}
	}
	return answer
}

func mulMtrMtr(A [m][n]float64, AT [n][m]float64) [n][n]float64 {
	var answer [n][n]float64
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < m; k++ {
				answer[i][j] += AT[i][k] * A[k][j]
			}
		}
	}
	return answer
}

func L(A_ [m][n]float64, u [n]float64, f_ [m]float64, h, sigma float64) float64 {
	var buff [3]float64
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			buff[0] += A_[i][j] * u[j]
		}
		buff[0] -= f_[i]
		buff[1] += buff[0] * buff[0]
	}
	buff[1] = math.Sqrt(buff[1])
	for j := 0; j < n; j++ {
		buff[2] += u[j] * u[j]
	}
	buff[2] = math.Sqrt(buff[2])
	return buff[1] + h*buff[2] + sigma
}

func F(A_ [m][n]float64, u [n]float64, f_ [m]float64, h float64) [n]float64 {
	var mass [n]float64
	for i := 0; i < n; i++ {
		mass[i] = diff(h, i, A_, u, f_)
	}
	return mass
}

func diff(h float64, k int, A_ [m][n]float64, u [n]float64, f_ [m]float64) float64 {
	var buff [5]float64
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			buff[0] += A_[i][j] * u[j]
		}
		buff[0] -= f_[i]
		buff[1] += A_[i][k] * buff[0]
		buff[2] += buff[0] * buff[0]
	}
	buff[2] = math.Sqrt(buff[2])
	buff[3] = h * u[k]
	for j := 0; j < n; j++ {
		buff[4] += u[j] * u[j]
	}
	buff[4] = math.Sqrt(buff[4])
	return buff[1]/buff[2] + buff[3]/buff[4]
}

func gradientDescent(A_ [m][n]float64, x0 [n]float64, f_ [m]float64, h, sigma, eps float64, M int) [n]float64 {
	var Fp [n]float64
	k := 0
	xPrev := x0
	x := xPrev
	tK := 3.0
	for {
		k += 1
		Fp = F(A_, xPrev, f_, h)
		for j := 0; j < n; j++ {
			x[j] = xPrev[j] - tK*Fp[j]
		}

		if L(A_, x, f_, h, sigma)-L(A_, xPrev, f_, h, sigma) < 0 {
			for j := 0; j < n; j++ {
				Fp[j] = x[j] - xPrev[j]
			}
			if normVecN(Fp) <= eps || k+1 > M {
				return x
			}
		} else {
			tK /= 2.
		}
		xPrev = x
	}
}

func TransposeMtr(A [m][n]float64) [n][m]float64 {
	var AT [n][m]float64
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			AT[j][i] = A[i][j]
		}
	}
	return AT
}

func LinearCG(A [n][n]float64, b [n]float64) [n]float64 {
	var xk, pk, curveX [n]float64
	rk := mulMtrVecNN(A, xk)
	for i := 0; i < n; i++ {
		rk[i] -= b[i]
		pk[i] -= rk[i]
	}
	rkNorm := normVecN(rk)
	numIter := 0
	for rkNorm > 1e-5 {
		apk := mulMtrVecNN(A, pk)
		rk2 := mulVecVecNN(rk, rk)
		alpha := rk2 / mulVecVecNN(pk, apk)
		for i := 0; i < n; i++ {
			xk[i] += alpha * pk[i]
			rk[i] += alpha * apk[i]
		}
		beta := mulVecVecNN(rk, rk) / rk2
		for i := 0; i < n; i++ {
			pk[i] = -rk[i] + beta*pk[i]
		}
		numIter += 1
		curveX = xk
		rkNorm = normVecN(rk)
	}
	return curveX
}

func rho(A_ [m][n]float64, u [n]float64, f_ [m]float64, lambdaDelta, sigma, h float64) float64 {
	var buff [5]float64

	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			buff[0] += A_[i][j] * u[j]
		}
		buff[0] -= f_[i]
		buff[1] += buff[0] * buff[0]
	}
	for j := 0; j < n; j++ {
		buff[2] += u[j] * u[j]
	}
	buff[2] = math.Sqrt(buff[2])
	return buff[1] - (lambdaDelta+sigma+h*buff[2])*(lambdaDelta+sigma+h*buff[2])
}

func main() {
	var A_, A, buffMtr [m][n]float64
	var u [n]float64
	var f, f_, buffVec [m]float64
	rand.Seed(1)
	for i := 0; i < n; i++ {
		u[i] = a + (b-a)*rand.Float64()
		for j := 0; j < n; j++ {
			A[i][j] = a + (b-a)*rand.Float64()
		}
	}
	f = mulMtrVecMN(A, u)
	for j := 0; j < n; j++ {
		A[n][j] = A[0][j] + A[n-1][j]
		A[n+1][j] = A[1][j] + A[n-2][j]
	}
	f[n] = f[0] + f[n-1]
	f[n+1] = f[1] + f[n-2]
	randA := -1.0
	randB := 1.0
	epsA := 0.001
	epsF := 0.001
	for i := 0; i < m; i++ {
		f_[i] = f[i] * (1 + (randA+(randB-randA)*rand.Float64())*epsF)
		buffVec[i] = f[i] - f_[i]
		for j := 0; j < n; j++ {
			A_[i][j] = A[i][j] * (1 + (randA+(randB-randA)*rand.Float64())*epsA)
			buffMtr[i][j] = A[i][j] - A_[i][j]
		}
	}
	h := normMtr(buffMtr)
	sigma := normVecM(buffVec)
	solution1 := gradientDescent(A_, u, f_, h, sigma, 1e-6, 500)
	lambdaDelta := L(A_, solution1, f_, h, sigma)
	AT := TransposeMtr(A_)
	alpha := 1.0
	val := 1.0
	rhoVec := make([]float64, 0)
	var answer [n]float64
	for val >= 1e-6 {
		ATA := mulMtrMtr(A_, AT)
		vec := mulMtrVecNM(AT, f_)
		for i := 0; i < n; i++ {
			ATA[i][i] += alpha
		}
		answer = LinearCG(ATA, vec)
		val = rho(A_, answer, f_, lambdaDelta, sigma, h)
		rhoVec = append(rhoVec, val)
		alpha -= 0.0001
	}
	createLineChart(rhoVec)
	fmt.Println(val, alpha)
	fmt.Println(answer)
	fmt.Println(u)
	fmt.Println(normVecN(answer))
	fmt.Println(normVecN(u))
}
