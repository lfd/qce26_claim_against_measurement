#let bra(f) = $lr(chevron.l #f|)$
#let ket(f) = $lr(|#f chevron.r)$

#let imag = $upright(i)$
#let euler = $upright(e)$

#let braket(..sink) = {
  let args = sink.pos() // array

  assert(1 <= args.len() and args.len() <= 3, message: "expecting 1, 2, or 3 args")

  let bra = args.at(0, default: none)
  let ket = args.at(-1, default: bra)

  if args.len() <= 2 {
    $ lr(chevron.l bra#h(0pt)mid(|)#h(0pt)ket chevron.r) $
  } else {
    let middle = args.at(1)
    $ lr(chevron.l bra#h(0pt)mid(|)#h(0pt)middle#h(0pt)mid(|)#h(0pt)ket chevron.r) $
  }
}

#let ketbra(..sink) = {
  let args = sink.pos() // array
  assert(args.len() == 1 or args.len() == 2, message: "expecting 1 or 2 args")

  let ket = args.at(0)
  let bra = args.at(1, default: ket)

  $ lr(|ket#h(0pt)mid(chevron.r#h(0pt)chevron.l)#h(0pt)bra|) $
}

#let matrixelement(n, M, m) = {
  $ lr(chevron.l #n#h(0pt)mid(|)#h(0pt)#M#h(0pt)mid(|)#h(0pt)#m chevron.r) $
}

#let x = $x$
#let y = $y$

#let NP = $upright(N P)$

#let tensor = $times.o$

#let H = $bold(H)$
#let G = $bold(G)$
#let c = $bold(c)$
#let e = $bold(e)$
#let s = $bold(s)$
#let m = $bold(m)$
#let x = $bold(x)$
#let y = $bold(y)$
#let v = $bold(v)$
#let B = $bold(B)$
#let b = $bold(b)$
#let l = $cal(l)$
#let A = $bold(A)$
#let a = $bold(a)$
#let zero = $bold(0)$

#let infobox(body) = block(
  fill: rgb("#a8abb3"),
  inset: 12pt,
  radius: 6pt,
  stroke: (paint: rgb("#75767c")),
)[
  #body
]