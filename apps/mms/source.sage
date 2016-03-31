x = var('x')
y = var('y')

phi = sin(x)

k = 0.75

s = - ( diff(1*phi, x) + diff(0*phi, y) ) + k*( diff(diff(phi, x), x) + diff(diff(phi, y), y) )

print s
