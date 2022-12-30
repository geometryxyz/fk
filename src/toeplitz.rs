pub struct Toeplitz<T: PartialEq> {
    pub(crate) elems: Vec<Vec<T>>,
}

impl<T: PartialEq> Toeplitz<T> {
    pub fn new(elems: Vec<Vec<T>>) -> Self {
        Self { elems }
    }

    pub fn is_toeplitz(elems: &Vec<Vec<T>>) {
        let n = elems.len();
        let m = elems[0].len();

        for j in 0..m {
            Self::check_diagonal(elems, 0, j);
        }

        for i in 1..n {
            Self::check_diagonal(elems, i, 0);
        }
    }

    fn check_diagonal(elems: &Vec<Vec<T>>, i_start: usize, j_start: usize) {
        let n = elems.len();
        let m = elems[0].len();

        for (i, j) in (i_start + 1..n).zip(j_start + 1..m) {
            if elems[i][j] != elems[i_start][j_start] {
                panic!("Not Toeplitz");
            }
        }
    }
}

#[cfg(test)]
mod toeplitz_test {
    use super::Toeplitz;

    #[test]
    fn test_check_correctness() {
        let elems = vec![
            vec![6, 7, 8, 9],
            vec![4, 6, 7, 8],
            vec![1, 4, 6, 7],
            vec![0, 1, 4, 6],
            vec![2, 0, 1, 4],
        ];

        Toeplitz::<usize>::is_toeplitz(&elems);
    }
}
