Lot na Księżyc
==============

Model matematyczny
------------------

### Równania rózniczkowe

$\dot{x}_1 = f_1(x(t),u(t),t) = x_3$

$\dot{x}_2 = f_2(x(t),u(t),t) = x_4$

$\dot{x}_3 = f_3(x(t),u(t),t) = x_1 + 2 \cdot x_4 -\frac{1 - \mu}{\gamma_E^3}
\cdot (x_1 + \mu) -\frac{\mu}{\gamma_M^3} \cdot (x_1 - 1 + \mu) + \frac{u_r
\cdot \cos(\phi)}{x_5}$

$\dot{x}_4 = f_4(x(t),u(t),t) = x_2 - 2 \cdot x_3 -\frac{1 - \mu}{\gamma_E^3}
\cdot x_2 -\frac{\mu_M}{\gamma_M^3} \cdot x_2 + \frac{u_r \cdot
\sin(\phi)}{x_5}$

$\dot{x}_5 = f_5(x(t),u(t),t) = -C_1 \cdot u_r$

Gdzie:

$\gamma_E = \sqrt{(x_1 + \mu)^2 + x_2^2}$

$\gamma_M = \sqrt{(x_1 - 1 + \mu)^2 + x_2^2}$

### Sterowanie

$u(t) = [u_r\ \ \phi]^T$

### Funkcja celu

$max\ Q_1$

$Q_1(u,x_0,T,x_T) = x_5(T) - K_1 \cdot T$

Równoważne:

$min\ Q_2$

$Q_2(u,x_0,T,x_T) = - x_5(T) + K_1 \cdot T$

### Warunki początkowe

$x_1(0) = r_E \cdot \cos(\theta_E) - \mu$

$x_2(0) = r_E \cdot \sin(\theta_E)$

$x_3(0) = \Delta V \cdot \sin(\alpha - \theta_E) - V_E \cdot \sin(\theta_E) +
x_2$

$x_4(0) = \Delta V \cdot \cos(\alpha - \theta_E) + V_E \cdot \cos(\theta_E) -
x_1$

$x_5(0)=m_0 \cdot e^{-C_2 \cdot \Delta V}$

Gdzie:

$r_E = \frac{a_E + LEO}{D}$

$V_E = \sqrt{\frac{1 - \mu}{r_E}}$

### Warunki końcowe

$r_M = \sqrt{(x_1(T) - 1 + \mu)^2+x_2(T)^2}$

$V_M = \sqrt{(x_3(T) - x_2(T))^2 + (x_4(T) + x_1(T))^2}$

$0 = (x_1(T) - 1 + \mu) \cdot (x_3(T) - x_2(T)) + x_2(T) \cdot (x_4(T) +
x_1(T))$

Gdzie:

$r_M = \frac{a_M + LMO}{D}$

$V_M = \sqrt{\frac{\mu}{r_M}}$

### Ograniczenia zmiennych stanu

$x_5(T) \ge m_r$

### Ograniczenia i warunki końcowe w postaci zmodyfikowanej

$k_1 = \frac{1}{4} \cdot ((x_1(T) - 1 + \mu)^2+x_2(T)^2 - r_M^2)^2$

$k_2 = \frac{1}{4} \cdot ((x_3(T) - x_2(T))^2 + (x_4(T) + x_1(T))^2 - V_M^2)^2$

$k_3 = \frac{1}{2} \cdot ((x_1(T) - 1 + \mu) \cdot (x_3(T) - x_2(T)) + x_2(T)
\cdot (x_4(T) + x_1(T)))^2$

$k_4 = \begin{cases} \frac{1}{2} \cdot (x_5(T) - m_r)^2, & x_5(T) < m_r \ 0, &
m_r \ge x_5(T) \end{cases}$

### Zewnętrzna funkcja kary

$\Phi(K) = K_2 \cdot (k_1 + k_2 + k_3 + k_4)$

### Zmodyfikowany wskaźnik jakości

$Q(u,x_0,T,x_T) = - x_5(T) + K_1 \cdot T + K_2 \cdot (k_1 + k_2 + k_3 + k_4)$

Zasada maksimum
---------------

### Hamiltonian

$H = \Psi_1 \cdot f_1 + \Psi_2 \cdot f_2 + \Psi_3 \cdot f_3 + \Psi_4 \cdot f_4 +
\Psi_5 \cdot f_5$

### Równania sprzężone

$\dot{\Psi}_1 = \Psi_3 \cdot \frac{\partial f_3}{\partial x_1} + \Psi_4 \cdot
\frac{\partial f_4}{\partial x_1}$

$\dot{\Psi}_2 = \Psi_3 \cdot \frac{\partial f_3}{\partial x_2} + \Psi_4 \cdot
\frac{\partial f_4}{\partial x_2}$

$\dot{\Psi}_3 = - \Psi_1 + 2 \cdot \Psi_4$

$\dot{\Psi}_4 = - \Psi_2 - 2 \cdot \Psi_3$

$\dot{\Psi}_5 = \Psi_3 \cdot \frac{u_r \cdot \cos(\phi)}{x_5^2} + \Psi_4 \cdot
\frac{u_r \cdot \sin(\phi)}{x_5^2}$

Gdzie:

$\frac{\partial f_3}{\partial x_1} = - 1 + \frac{(1 - \mu) \cdot (\gamma_E^2 - 3
\cdot (x_1 + \mu)^2)}{\gamma_E^5} + \frac{\mu \cdot (\gamma_M^2 - 3 \cdot (x_1 -
1 + \mu)^2)}{\gamma_M^5}$

$\frac{\partial f_4}{\partial x_1} = -\frac{3 \cdot (1 - \mu) \cdot (x_1 + \mu)
\cdot x_2}{\gamma_E^5} - \frac{3 \cdot \mu \cdot (x_1 - 1 + \mu) \cdot
x_2}{\gamma_M^5}$

$\frac{\partial f_3}{\partial x_2} = = -\frac{3 \cdot (1 - \mu) \cdot (x_1 +
\mu) \cdot x_2}{\gamma_E^5} - \frac{3 \cdot \mu \cdot (x_1 - 1 + \mu) \cdot
x_2}{\gamma_M^5}$

$\frac{\partial f_4}{\partial x_2} = - 1 + \frac{(1 - \mu) \cdot (\gamma_E^2 - 3
\cdot x_2^2)}{\gamma_E^5} + \frac{\mu \cdot (\gamma_M^2 - 3 \cdot
x_2^2)}{\gamma_M^5}$

### Warunki końcowe dla równań sprzężonych

$\Psi_1(T) = K_2 \cdot (\beta_1 \cdot (x_1 - 1 + \mu) + \beta_2 \cdot (x_4 +
x_1) + \beta_3 \cdot x_3(T))$

$\Psi_2(T) = K_2 \cdot (\beta_1 \cdot x_2 + \beta_2 \cdot (x_2(T) - x_3(T)) +
\beta_3 \cdot (1 - \mu + x_4(T)))$

$\Psi_3(T) = K_2 \cdot (\beta_2 \cdot (x_3(T) - x_2(T)) + \beta_3 \cdot (x_1 -1
+ \mu))$

$\Psi_4(T) = K_2 \cdot (\beta_2 \cdot (x_4(T) + x_1(T)) + \beta_3 \cdot x_2(T))$

$\Psi_5(T) = 1 - K_2 \cdot \begin{cases} x_5(T) - m_r, & x_5(T) \le m_r \ 0, &
m_r < x_5(T) \end{cases}$

Gdzie:

$\beta_1 = r_M^2 - (x_1(T) - 1 + \mu)^2 - x_2(T)^2$

$\beta_2 = V_M^2 - (x_3(T) - x_2(T))^2 - (x_4(T)+ x_1(T))^2$

$\beta_3 = - (x_1(T) - 1 + \mu) \cdot (x_3(T) - x_2(T)) - x_2(T) \cdot (x_4(T) +
x_1(T))$

### Źródła

<http://www.hindawi.com/journals/mpe/2012/971983/>

Welcome to Texts!
=================

Texts is “rich editor for plain text”—you can write documents in visual
interface and they will be saved in plain text format with Markdown markup. This
document shows some styles available in Texts. Please feel free to edit it—a
fresh copy is always available via Help menu.

Inline Formatting
-----------------

To apply *emphasis* or **strong emphasis** press Control+I or Control+B and
type. Press again to turn it off.

Tables
------

+---+---+---+
| 1 | 2 | 3 |
+---+---+---+
| 4 | 5 | 6 |
+---+---+---+
| 7 | 8 | 9 |
+---+---+---+

Footnotes
---------

Notes can be placed anywhere[^1] in the source file (press Control+R to insert
one). When document is published to paper-oriented format (PDF or Word)
footnotes are displayed at the bottom of the corresponding page. E-book readers
usually display notes in popups.

[^1]: Yes, right here.

Hyperlinks
----------

Press Control+K to insert a hyperlink. It can be either literal URL
(<http://www.google.com/>) or have some [text][2]. Click URL with Control key
pressed to open it in web browser.

[2]: <http://www.texts.io/>

Math
----

Formulas can be placed inline like $E=mc^2$ or in a separate paragraph, like the
following one. Standard LaTeX syntax is supported.

$$
\frac{n!}{k!(n-k)!} = \binom{n}{k}
$$

Code
----

Inline `code` gets monospaced font.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Verbatim blocks use monospaced font as well and preserve line
breaks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Press “Enter” key inside code block to insert a line break or “Shift+Enter” to
end code block.

Lists
-----

-   First bulleted item.

-   Second bulleted item.

    -   Subitem.

        -   Subsubitem.

            -   And some more.

1.  First numbered item.

2.  Second numbered item.

    1.  Subitem.

        1.  Subsubitem.

Lists can be entered via autoformatting: type minus and space for bulleted item
or “1”, point and space for numbered item. Press “Tab” or “Shift+Tab” to
indent/unindent items.

Happy writing!
==============

Got a question? Visit <http://www.texts.io/support/>.
