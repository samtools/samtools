/*  misc/HmmGlocal.java.

    Copyright (C) 2010 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

import java.io.*;
import java.lang.*;

public class HmmGlocal
{
    private double[] qual2prob;
    private double cd, ce; // gap open probility [1e-3], gap extension probability [0.1]
    private int cb; // band width [7]

    public HmmGlocal(final double d, final double e, final int b) {
        cd = d; ce = e; cb = b;
        qual2prob = new double[256];
        for (int i = 0; i < 256; ++i)
            qual2prob[i] = Math.pow(10, -i/10.);
    }
    private static int set_u(final int b, final int i, final int k) {
        int x = i - b;
        x = x > 0? x : 0;
        return (k + 1 - x) * 3;
    }
    public int hmm_glocal(final byte[] _ref, final byte[] _query, final byte[] _iqual, int[] state, byte[] q) {
        int i, k;
        /*** initialization ***/
        // change coordinates
        int l_ref = _ref.length;
        byte[] ref = new byte[l_ref+1];
        for (i = 0; i < l_ref; ++i) ref[i+1] = _ref[i]; // FIXME: this is silly...
        int l_query = _query.length;
        byte[] query = new byte[l_query+1];
        double[] qual = new double[l_query+1];
        for (i = 0; i < l_query; ++i) {
            query[i+1] = _query[i];
            qual[i+1] = qual2prob[_iqual[i]];
        }
        // set band width
        int bw2, bw = l_ref > l_query? l_ref : l_query;
        if (bw > cb) bw = cb;
        if (bw < Math.abs(l_ref - l_query)) bw = Math.abs(l_ref - l_query);
        bw2 = bw * 2 + 1;
        // allocate the forward and backward matrices f[][] and b[][] and the scaling array s[]
        double[][] f = new double[l_query+1][bw2*3 + 6];
        double[][] b = new double[l_query+1][bw2*3 + 6];
        double[] s = new double[l_query+2];
        // initialize transition probabilities
        double sM, sI, bM, bI;
        sM = sI = 1. / (2 * l_query + 2);
        bM = (1 - cd) / l_query; bI = cd / l_query; // (bM+bI)*l_query==1
        double[] m = new double[9];
        m[0*3+0] = (1 - cd - cd) * (1 - sM); m[0*3+1] = m[0*3+2] = cd * (1 - sM);
        m[1*3+0] = (1 - ce) * (1 - sI); m[1*3+1] = ce * (1 - sI); m[1*3+2] = 0.;
        m[2*3+0] = 1 - ce; m[2*3+1] = 0.; m[2*3+2] = ce;
        /*** forward ***/
        // f[0]
        f[0][set_u(bw, 0, 0)] = s[0] = 1.;
        { // f[1]
            double[] fi = f[1];
            double sum;
            int beg = 1, end = l_ref < bw + 1? l_ref : bw + 1, _beg, _end;
            for (k = beg, sum = 0.; k <= end; ++k) {
                int u;
                double e = (ref[k] > 3 || query[1] > 3)? 1. : ref[k] == query[1]? 1. - qual[1] : qual[1] / 3.;
                u = set_u(bw, 1, k);
                fi[u+0] = e * bM; fi[u+1] = .25 * bI;
                sum += fi[u] + fi[u+1];
            }
            // rescale
            s[1] = sum;
            _beg = set_u(bw, 1, beg); _end = set_u(bw, 1, end); _end += 2;
            for (k = _beg; k <= _end; ++k) fi[k] /= sum;
        }
        // f[2..l_query]
        for (i = 2; i <= l_query; ++i) {
            double[] fi = f[i], fi1 = f[i-1];
            double sum, qli = qual[i];
            int beg = 1, end = l_ref, x, _beg, _end;
            byte qyi = query[i];
            x = i - bw; beg = beg > x? beg : x; // band start
            x = i + bw; end = end < x? end : x; // band end
            for (k = beg, sum = 0.; k <= end; ++k) {
                int u, v11, v01, v10;
                double e;
                e = (ref[k] > 3 || qyi > 3)? 1. : ref[k] == qyi? 1. - qli : qli / 3.;
                u = set_u(bw, i, k); v11 = set_u(bw, i-1, k-1); v10 = set_u(bw, i-1, k); v01 = set_u(bw, i, k-1);
                fi[u+0] = e * (m[0] * fi1[v11+0] + m[3] * fi1[v11+1] + m[6] * fi1[v11+2]);
                fi[u+1] = .25 * (m[1] * fi1[v10+0] + m[4] * fi1[v10+1]);
                fi[u+2] = m[2] * fi[v01+0] + m[8] * fi[v01+2];
                sum += fi[u] + fi[u+1] + fi[u+2];
                //System.out.println("("+i+","+k+";"+u+"): "+fi[u]+","+fi[u+1]+","+fi[u+2]);
            }
            // rescale
            s[i] = sum;
            _beg = set_u(bw, i, beg); _end = set_u(bw, i, end); _end += 2;
            for (k = _beg, sum = 1./sum; k <= _end; ++k) fi[k] *= sum;
        }
        { // f[l_query+1]
            double sum;
            for (k = 1, sum = 0.; k <= l_ref; ++k) {
                int u = set_u(bw, l_query, k);
                if (u < 3 || u >= bw2*3+3) continue;
                sum += f[l_query][u+0] * sM + f[l_query][u+1] * sI;
            }
            s[l_query+1] = sum; // the last scaling factor
        }
        /*** backward ***/
        // b[l_query] (b[l_query+1][0]=1 and thus \tilde{b}[][]=1/s[l_query+1]; this is where s[l_query+1] comes from)
        for (k = 1; k <= l_ref; ++k) {
            int u = set_u(bw, l_query, k);
            double[] bi = b[l_query];
            if (u < 3 || u >= bw2*3+3) continue;
            bi[u+0] = sM / s[l_query] / s[l_query+1]; bi[u+1] = sI / s[l_query] / s[l_query+1];
        }
        // b[l_query-1..1]
        for (i = l_query - 1; i >= 1; --i) {
            int beg = 1, end = l_ref, x, _beg, _end;
            double[] bi = b[i], bi1 = b[i+1];
            double y = (i > 1)? 1. : 0., qli1 = qual[i+1];
            byte qyi1 = query[i+1];
            x = i - bw; beg = beg > x? beg : x;
            x = i + bw; end = end < x? end : x;
            for (k = end; k >= beg; --k) {
                int u, v11, v01, v10;
                double e;
                u = set_u(bw, i, k); v11 = set_u(bw, i+1, k+1); v10 = set_u(bw, i+1, k); v01 = set_u(bw, i, k+1);
                e = (k >= l_ref? 0 : (ref[k+1] > 3 || qyi1 > 3)? 1. : ref[k+1] == qyi1? 1. - qli1 : qli1 / 3.) * bi1[v11];
                bi[u+0] = e * m[0] + .25 * m[1] * bi1[v10+1] + m[2] * bi[v01+2]; // bi1[v11] has been foled into e.
                bi[u+1] = e * m[3] + .25 * m[4] * bi1[v10+1];
                bi[u+2] = (e * m[6] + m[8] * bi[v01+2]) * y;
            }
            // rescale
            _beg = set_u(bw, i, beg); _end = set_u(bw, i, end); _end += 2;
            for (k = _beg, y = 1./s[i]; k <= _end; ++k) bi[k] *= y;
        }
        double pb;
        { // b[0]
            int beg = 1, end = l_ref < bw + 1? l_ref : bw + 1;
            double sum = 0.;
            for (k = end; k >= beg; --k) {
                int u = set_u(bw, 1, k);
                double e = (ref[k] > 3 || query[1] > 3)? 1. : ref[k] == query[1]? 1. - qual[1] : qual[1] / 3.;
                if (u < 3 || u >= bw2*3+3) continue;
                sum += e * b[1][u+0] * bM + .25 * b[1][u+1] * bI;
            }
            pb = b[0][set_u(bw, 0, 0)] = sum / s[0]; // if everything works as is expected, pb == 1.0
        }
        int is_diff = Math.abs(pb - 1.) > 1e-7? 1 : 0;
        /*** MAP ***/
        for (i = 1; i <= l_query; ++i) {
            double sum = 0., max = 0.;
            double[] fi = f[i], bi = b[i];
            int beg = 1, end = l_ref, x, max_k = -1;
            x = i - bw; beg = beg > x? beg : x;
            x = i + bw; end = end < x? end : x;
            for (k = beg; k <= end; ++k) {
                int u = set_u(bw, i, k);
                double z;
                sum += (z = fi[u+0] * bi[u+0]); if (z > max) { max = z; max_k = (k-1)<<2 | 0; }
                sum += (z = fi[u+1] * bi[u+1]); if (z > max) { max = z; max_k = (k-1)<<2 | 1; }
            }
            max /= sum; sum *= s[i]; // if everything works as is expected, sum == 1.0
            if (state != null) state[i-1] = max_k;
            if (q != null) {
                k = (int)(-4.343 * Math.log(1. - max) + .499);
                q[i-1] = (byte)(k > 100? 99 : k);
            }
            //System.out.println("("+pb+","+sum+")"+" ("+(i-1)+","+(max_k>>2)+","+(max_k&3)+","+max+")");
        }
        return 0;
    }

    public static void main(String[] args) {
        byte[] ref = {'\0', '\1', '\3', '\3', '\1'};
        byte[] query = {'\0', '\3', '\3', '\1'};
        byte[] qual = new byte[4];
        qual[0] = qual[1] = qual[2] = qual[3] = (byte)20;
        HmmGlocal hg = new HmmGlocal(1e-3, 0.1, 7);
        hg.hmm_glocal(ref, query, qual, null, null);
    }
}
