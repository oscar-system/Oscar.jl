import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: "A VitePress Site",
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }]],
  ignoreDeadLinks: true,
  appearance: 'dark',

  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },
  themeConfig: {
    outline: 'deep',
    logo: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Getting Started', link: '/start/' },
      { text: 'Manual', link: '/manual/' },
      { text: 'Tutorials', link: '/tutorials/' },
      { text: 'How-to Guides', link: '/howto/' },
      {
        text: 'Versions', items: [
          { text: 'Stable', link: 'https://thofma.com/Hecke.jl/stable/' },
          { text: 'Dev', link: 'https://thofma.com/Hecke.jl/dev/' }
        ]
      }
    ],
    sidebar: {
      '/manual/': [
        { text: 'Manual', 
          items: [
            { text: 'Introduction', link: '/manual/'},
            { text: 'Number fields',
              collapsed: true,
              items: [
                { text: 'Introduction', link: 'manual/number_fields/intro'},
                { text: 'Number field operations', link: 'manual/number_fields/fields'},
                { text: 'Elements', link: 'manual/number_fields/elements'},
                { text: 'Complex embeddings', link: 'manual/number_fields/complex_embeddings'},
                { text: 'Class field theory', link: 'manual/number_fields/class_fields.md'}, 
                { text: 'Internals', link: 'manual/number_fields/internal'}
              ]
            },
            {
              text: 'Orders in number fields',
              collapsed: true,
              items: [
                { text: 'Introduction', link: 'manual/orders/introduction'},
                { text: 'Basics', link: 'manual/orders/orders'},
                { text: 'Elements', link: 'manual/orders/elements'},
                { text: 'Ideals', link: 'manual/orders/ideals'},
                { text: 'Fractional ideals', link: 'manual/orders/frac_ideals'},
              ]
            },
            {
              text: 'Quadratic and hermitian forms',
              collapsed: true,
              items: [
                { text: 'Introduction', link: 'manual//quad_forms/introduction'},
                { text: 'Basics', link: 'manual//quad_forms/basics'},
                { text: 'Lattices', link: 'manual/quad_forms/lattices'},
                { text: 'Integer lattices', link: 'manual/quad_forms/integer_lattices'},
                { text: 'Genera for hermitian lattices', link: 'manual/quad_forms/genusherm'},
                { text: 'Genera for integer lattices', link: 'manual/quad_forms/Zgenera'},
                { text: 'Discriminant groups', link: 'manual/quad_forms/discriminant_group'},
              ]
            },
            { text: 'Algebras',
              collapsed: true,
              items: [
                { text: 'Introduction', link: 'manual/algebras/intro'},
                { text: 'Basics', link: 'manual/algebras/basics'},
                { text: 'Structure constant algebras', link: 'manual/algebras/structureconstant'},
              ]
            },
            {
              text: 'Elliptic curves',
              collapsed: true,
              items: [
                { text: 'Introduction', link: 'manual/elliptic_curves/intro'},
                { text: 'Basics', link: 'manual/elliptic_curves/basics'},
                { text: 'Finite fields', link: 'manual/elliptic_curves/finite_fields'},
                { text: 'Number fields', link: 'manual/elliptic_curves/number_fields'},
              ]
            },
            {
              text: 'Abelian groups',
              collapsed: true,
              items: [
                { text: 'Introduction', link: 'manual/abelian/introduction'},
                { text: 'Elements', link: 'manual/abelian/elements'},
                { text: 'Operations', link: 'manual/abelian/structural'},
                { text: 'Morphisms', link: 'manual/abelian/maps'},
              ]
            },
            {
              text: 'Miscellaneous',
              collapsed: true,
              items: [
                { text: 'Factored elements', link: 'manual/misc/FacElem'},
                { text: 'Sparse linear algebra', link: 'manual/misc/sparse'},
                { text: 'Conjugacy of integer matrices', link: 'manual/misc/conjugacy'},
                { text: 'Multisets', link: 'manual/misc/mset'},
                { text: 'Pseudo-matrices', link: 'manual/misc/pmat'},
              ]
            },
            {
              text: 'Developer',
              collapsed: true,
              items: [
                { text: 'Macros', link: 'manual/developer/macros'},
                { text: 'Tests', link: 'manual/developer/test'},
              ]
            }
          ]
        }
      ],
      '/tutorials/': [
        { text: 'Tutorials', 
          items: [
            { text: 'Add me', link: '/tutorials/'},
          ]
        }
      ],
      '/howto/': [
        { text: 'How-to Guides', 
          items: [
            { text: 'Reduction of polynomials', link: '/howto/reduction'},
          ]
        }
      ],

    },
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
